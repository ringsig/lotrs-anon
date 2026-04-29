"""
lotrs.py -- LoTRS: Practical Post-Quantum Threshold Ring Signatures
             from Lattices.

Reference implementation of the scheme described in Figs. 4-6 of the
LoTRS paper.  Every source of randomness is derived deterministically
from caller-supplied seeds via SHAKE-256, enabling known-answer testing.

Notation
--------
- R_q = Z_q[X] / (X^d + 1)         main ring
- R_qhat = Z_{q_hat}[X] / (X^d+1)  binary-proof ring
- PK  = N columns x T rows  of public keys  t_{u,i} in R_q^k
- ell = hidden column index  in [0, N-1]
"""

import math

from ring import Ring
from sample import (build_cdt, make_xof, derive_subseed,
                    build_cdt_sampler, build_facct_sampler,
                    xof_sample_uniform, xof_sample_gaussian,
                    xof_sample_gaussian_facct,
                    xof_sample_short, xof_sample_challenge,
                    xof_sample_bounded, rej, rej_op, _flat)
from params import LoTRSParams


class LoTRS:
    """LoTRS structured threshold ring signature scheme.

    This implementation covers the kappa = 1 case only (the setting
    used in the concrete parameter set of the paper).  Attempting to
    instantiate with kappa > 1 will raise.
    """

    def __init__(self, par: LoTRSParams):
        par.check()
        if par.kappa != 1:
            raise NotImplementedError(
                "only kappa = 1 is implemented; "
                f"got kappa = {par.kappa}")
        self.par = par
        self.Rq = Ring(par.q, par.d)
        self.Rqh = Ring(par.q_hat, par.d)

        # Gaussian-sampler layout — no dispatch on sigma here; the
        # mask backend is declared in `par.mask_sampler`.  sigma_a /
        # sigma_b are binary-proof widths and always use a CDT.
        self.cdt_a = build_cdt(par.sigma_a, par.lam)
        self.cdt_b = build_cdt(par.sigma_b, par.lam)

        if par.mask_sampler == "cdt":
            self.gauss_0 = build_cdt_sampler(par.sigma_0, par.lam)
            self.gauss_0p = build_cdt_sampler(par.sigma_0_prime, par.lam)
        elif par.mask_sampler == "facct":
            self.gauss_0 = build_facct_sampler(par.sigma_0, par.lam)
            self.gauss_0p = build_facct_sampler(par.sigma_0_prime, par.lam)
        else:
            raise ValueError(
                f"unsupported par.mask_sampler = {par.mask_sampler!r}")

    # ==================================================================
    #  Fig. 4  --  Setup
    # ==================================================================

    def setup(self, seed):
        """
        Setup(1^lambda) -> pp.           [Fig. 4, Setup]

        pp is the 32-byte public seed  rho_pp  from which all public
        matrices (A, G) are deterministically expanded.
        """
        assert len(seed) == 32
        return seed                       # pp := rho_pp

    # ==================================================================
    #  Fig. 4  --  KGen
    # ==================================================================

    def keygen(self, pp, seed):
        """
        KGen(pp) -> (sk, pk).            [Fig. 4, KGen]

        Args
        ----
        pp   : public parameters  (32-byte seed)
        seed : 32-byte per-user key-generation seed

        Returns
        -------
        sk : list of  (l+k)  ring elements   (secret key  s)
        pk : list of  k      ring elements   (public key  t = A_bar * s)
        """
        par = self.par
        Rq = self.Rq

        # 1: A := Expand(rho_pp, A; R_q^{k x l})
        A = self._expand_A(pp)
        # 2: A_bar := [A | I_{k}]
        A_bar = self._augment_I(A, par.k)

        # 3: s <- S_eta^{l+k}
        xof = make_xof(seed, b"keygen")
        s = [xof_sample_short(xof, par.eta, par.d)
             for _ in range(par.l + par.k)]
        s = [Rq.from_centered(si) for si in s]

        # 4: t := A_bar * s
        t = Rq.mat_vec(A_bar, s)

        return s, t

    # ==================================================================
    #  Fig. 4  --  KAgg
    # ==================================================================

    def kagg(self, pk_table, pk_hash=None):
        """
        KAgg(PK) -> [t_tilde_0 .. t_tilde_{N-1}].   [Fig. 4, KAgg]

        pk_table[col_i][row_u] = t_{u,i}  (list of k ring elements)
        """
        par, Rq = self.par, self.Rq
        N, T, k = par.N, par.T, par.k

        if pk_hash is None:
            pk_hash = self._pk_hash(pk_table)
        alphas = [self._hash_agg(pk_hash, u) for u in range(T)]
        agg_keys = []
        for i in range(N):
            acc = Rq.vec_zero(k)
            for u in range(T):
                scaled = Rq.vec_scale(alphas[u], pk_table[i][u])
                acc = Rq.vec_add(acc, scaled)
            agg_keys.append(acc)
        return agg_keys

    # ==================================================================
    #  Fig. 4  --  Sign_1
    # ==================================================================

    def sign1(self, pp, sk_u, row_u, ell, mu, pk_table, rho, attempt,
              pk_hash=None):
        """
        Sign_1 for signer u.             [Fig. 4, Sign_1]

        Returns (state_u, Com_u)  where
        - state_u  carries all values needed by sign2
        - Com_u = (w_{u,j})_{j in [0, kappa-1]}  -- list of kappa
          vectors of k ring elements
        """
        par, Rq = self.par, self.Rq
        kappa, d = par.kappa, par.d
        if pk_hash is None:
            pk_hash = self._pk_hash(pk_table)

        # -- line 2: alpha_u
        alpha_u = self._hash_agg(pk_hash, row_u)

        # -- line 3: expand (a_{j,i})  from rho
        a_coeffs = self._expand_a_coeffs(rho, attempt)

        # -- selector coefficients  p_{i,j}  (Remark 2, p. 15)
        p_coeffs = self._compute_p_coeffs(a_coeffs, ell)

        # -- line 7: B := H_com(PK, mu)
        B = self._hash_com(pk_hash, mu)
        B_bar = self._augment_I(B, par.k)

        # -- line 1: A (from pp)
        A = self._expand_A(pp)
        A_bar = self._augment_I(A, par.k)

        # -- lines 8-10: sample masking vectors  y_{u,j}, r_{u,j}
        # Explicit dispatch on the sampler kind — keeps CDT and FACCT
        # call sites side-by-side instead of hiding them behind a
        # dispatcher helper.
        y_list, r_list = [], []
        for j in range(kappa):
            xof_y = make_xof(rho, b"y", row_u, j, attempt)
            if self.gauss_0.kind == "cdt":
                y_j = [xof_sample_gaussian(xof_y, self.gauss_0.data,
                                           par.lam, d)
                       for _ in range(par.l + par.k)]
            else:   # "facct"
                y_j = [xof_sample_gaussian_facct(xof_y, self.gauss_0.data, d)
                       for _ in range(par.l + par.k)]
            y_j = [Rq.from_centered(c) for c in y_j]
            y_list.append(y_j)

            xof_r = make_xof(rho, b"r", row_u, j, attempt)
            if self.gauss_0p.kind == "cdt":
                r_j = [xof_sample_gaussian(xof_r, self.gauss_0p.data,
                                           par.lam, d)
                       for _ in range(par.l_prime + par.k)]
            else:   # "facct"
                r_j = [xof_sample_gaussian_facct(xof_r, self.gauss_0p.data, d)
                       for _ in range(par.l_prime + par.k)]
            r_j = [Rq.from_centered(c) for c in r_j]
            r_list.append(r_j)

        # -- line 12: commitments  w_{u,j}
        #   w_{u,j} = A_bar y_{u,j} + B_bar r_{u,j}
        #             + alpha_u * sum_{i=0}^{N-1} p_{i,j} t_{u,i}
        coms = []
        for j in range(kappa):
            term_A = Rq.mat_vec(A_bar, y_list[j])
            term_B = Rq.mat_vec(B_bar, r_list[j])
            # public-key contribution
            pk_term = Rq.vec_zero(par.k)
            for i in range(par.N):
                # p_{i,j} is a ring element; t_{u,i} = pk_table[i][row_u]
                p_ij = p_coeffs[i][j]
                contrib = Rq.vec_scale(p_ij, pk_table[i][row_u])
                pk_term = Rq.vec_add(pk_term, contrib)
            pk_term = Rq.vec_scale(alpha_u, pk_term)

            w_uj = Rq.vec_add(Rq.vec_add(term_A, term_B), pk_term)
            coms.append(w_uj)

        state = dict(
            ell=ell, mu=mu, pp=pp, sk_u=sk_u, alpha_u=alpha_u,
            y_list=y_list, r_list=r_list, pk_table=pk_table, rho=rho,
            attempt=attempt, row_u=row_u, pk_hash=pk_hash,
        )
        return state, coms

    # ==================================================================
    #  Fig. 5  --  Sign_2
    # ==================================================================

    def sign2(self, state, all_coms, pk_table):
        """
        Sign_2 for signer u.             [Fig. 5, Sign_2 + Sign_bin]

        all_coms[u'][j] = w_{u',j}  -- round-1 commitments from every
        signer.

        Returns sigma_u = (pi, z_u, r_u)  or None on restart.
        """
        par, Rq = self.par, self.Rq
        kappa = par.kappa
        ell = state["ell"]
        mu = state["mu"]
        pp = state["pp"]
        rho = state["rho"]
        attempt = state["attempt"]
        row_u = state["row_u"]

        # -- lines 4-7: aggregate commitments  w_tilde_j
        w_tilde = []
        for j in range(kappa):
            acc = Rq.vec_zero(par.k)
            for u_prime in range(par.T):
                acc = Rq.vec_add(acc, all_coms[u_prime][j])
            w_tilde.append(acc)

        # -- line 5: decompose  w_tilde_j = 2^{K_w} w_tilde_j^(1)
        #                                  + w_tilde_j^(0)
        w_tilde_hi = []
        w_tilde_lo = []
        for j in range(kappa):
            hi_j, lo_j = [], []
            for comp in w_tilde[j]:
                h, lo = Rq.centered_decompose(comp, par.K_w)
                hi_j.append(h)
                lo_j.append(lo)
            w_tilde_hi.append(hi_j)
            w_tilde_lo.append(lo_j)

        # -- line 7: run binary selection proof  Sign_bin
        pi = self._sign_bin(pp, ell, mu, w_tilde_hi, pk_table,
                            rho, attempt, state["pk_hash"])
        if pi is None:
            return None                   # Sign_bin triggered Restart

        # extract Fiat-Shamir challenge  x  from pi
        x = pi["x"]

        # -- lines 9-11: w_tilde_0 stability check  [new paper]
        #   M_w = sum_{j>=1} ||x^j||_1 * 2^{K_{w,j}-1}
        #   For kappa=1, M_w = 0.
        M_w = 0
        if kappa > 1:
            x_pow_stab = [Rq.one()]
            for _ in range(kappa - 1):
                x_pow_stab.append(Rq.mul(x_pow_stab[-1], x))
            for j in range(1, kappa):
                M_w += Rq.l1_norm(x_pow_stab[j]) * (1 << (par.K_w - 1))
        threshold_w = (1 << (par.K_w - 1)) - M_w
        if Rq.vec_inf_norm(w_tilde_lo[0]) > threshold_w:
            return None                   # Restart

        # -- lines 9-11: main response  z_u
        #   z_u = x^kappa * alpha_u * s_u  -  sum_{j=0}^{kappa-1} x^j y_{u,j}
        alpha_u = state["alpha_u"]
        sk_u = state["sk_u"]
        y_list = state["y_list"]
        r_list = state["r_list"]

        x_pow = [Rq.one()]               # x^0, x^1, ..., x^kappa
        for _ in range(kappa):
            x_pow.append(Rq.mul(x_pow[-1], x))

        # x^kappa * alpha_u * s_u
        coeff = Rq.mul(x_pow[kappa], alpha_u)
        z_u = Rq.vec_scale(coeff, sk_u)

        for j in range(kappa):
            term = Rq.vec_scale(x_pow[j], y_list[j])
            z_u = Rq.vec_sub(z_u, term)

        # rejection sampling on z_u
        z_u_signed = [Rq.centered(c) for c in z_u]
        # shift  v = x^kappa alpha_u s_u - sum_{j>=1} x^j y_{u,j}
        #   for kappa=1:  v = x alpha_u s_u
        shift = Rq.vec_scale(coeff, sk_u)
        for j in range(1, kappa):
            shift = Rq.vec_sub(shift, Rq.vec_scale(x_pow[j], y_list[j]))
        shift_signed = [Rq.centered(c) for c in shift]

        xof_rej = make_xof(rho, b"rej_z", row_u, attempt)
        if not rej(xof_rej,
                   _flat(z_u_signed), _flat(shift_signed),
                   par.phi, par.B_0):
            return None                   # reject -> Restart

        # -- line 12: auxiliary response  r_u = - sum x^j r_{u,j}
        r_u = Rq.vec_zero(par.l_prime + par.k)
        for j in range(kappa):
            term = Rq.vec_scale(x_pow[j], r_list[j])
            r_u = Rq.vec_sub(r_u, term)

        return dict(pi=pi, z_u=z_u, r_u=r_u)

    # ==================================================================
    #  Fig. 4  --  SAgg
    # ==================================================================

    def sagg(self, sigmas):
        """
        SAgg({sigma_0, ..., sigma_{T-1}}) -> sigma_tilde.
                                              [Fig. 4, SAgg]
        """
        par, Rq = self.par, self.Rq
        T = par.T
        l, l_p, k = par.l, par.l_prime, par.k

        # -- input shape validation --
        # sigmas must be a sequence of exactly T dicts, each containing
        # the keys pi / z_u / r_u.  Malformed / dishonest input raises
        # ValueError rather than crashing with IndexError / KeyError.
        try:
            n_sigmas = len(sigmas)
        except TypeError:
            raise ValueError(
                f"SAgg: expected a sequence of signer transcripts, "
                f"got {type(sigmas).__name__}")
        if n_sigmas != T:
            raise ValueError(
                f"SAgg: expected {T} signer transcripts, got {n_sigmas}")

        required = ("pi", "z_u", "r_u")
        for u, sig_u in enumerate(sigmas):
            if not isinstance(sig_u, dict):
                raise ValueError(
                    f"SAgg: signer {u} transcript is not a dict "
                    f"(got {type(sig_u).__name__})")
            for key in required:
                if key not in sig_u:
                    raise ValueError(
                        f"SAgg: signer {u} transcript missing key {key!r}")

        # -- consistency check on pi (lines 2-6 of SAgg) --
        pi_required = ("x", "B_bin_hi", "w_tilde_hi", "f1", "z_b")
        pi0 = sigmas[0]["pi"]
        for key in pi_required:
            if key not in pi0:
                raise ValueError(
                    f"SAgg: signer 0 pi missing key {key!r}")

        for u in range(1, T):
            pi_u = sigmas[u]["pi"]
            for field in pi_required:
                if field not in pi_u:
                    raise ValueError(
                        f"SAgg: signer {u} pi missing key {field!r}")
                if pi_u[field] != pi0[field]:
                    raise ValueError(
                        f"SAgg: signer {u} disagrees with signer 0 "
                        f"on pi[{field!r}]")

        # split z_u  into  z'_u in R^l  and  z''_u in R^k
        # split r_u  into  r'_u in R^{l'} and  r''_u in R^k
        z_tilde = Rq.vec_zero(l)
        r_tilde = Rq.vec_zero(l_p)
        e_tilde = Rq.vec_zero(k)

        for u in range(T):
            z_u = sigmas[u]["z_u"]         # length l+k
            r_u = sigmas[u]["r_u"]         # length l'+k

            z_prime = z_u[:l]
            z_dprime = z_u[l:]
            r_prime = r_u[:l_p]
            r_dprime = r_u[l_p:]

            z_tilde = Rq.vec_add(z_tilde, z_prime)
            r_tilde = Rq.vec_add(r_tilde, r_prime)
            e_tilde = Rq.vec_add(e_tilde,
                                 Rq.vec_add(z_dprime, r_dprime))

        return dict(pi=pi0,
                    z_tilde=z_tilde,
                    r_tilde=r_tilde,
                    e_tilde=e_tilde)

    # ==================================================================
    #  Fig. 6  --  Vf  (verification)
    # ==================================================================

    def verify(self, pp, mu, sigma, pk_table):
        """
        Vf(pp, mu, sigma_tilde, PK) -> bool.   [Fig. 6]
        """
        par = self.par
        Rq, Rqh = self.Rq, self.Rqh
        kappa, d = par.kappa, par.d

        pi = sigma["pi"]
        z_tilde = sigma["z_tilde"]
        r_tilde = sigma["r_tilde"]
        e_tilde = sigma["e_tilde"]
        pk_hash = self._pk_hash(pk_table)

        # -- expand G and extract pi fields
        B_bin_hi = pi["B_bin_hi"]
        w_tilde_hi_rest = pi["w_tilde_hi"]   # j >= 1 only (kappa-1 entries)
        x = pi["x"]
        f1 = pi["f1"]
        z_b = pi["z_b"]

        # -- line 3: G
        G = self._expand_G(pp)

        # -- line 4: norm bounds on f1, f0 (reconstructed), z_b.
        # f0 is reconstructed from the challenge and f1 below; the
        # paper's Vf (Fig. 6) checks both ||f_1||_inf <= B_f1 and
        # ||f_0||_inf <= B_f0.
        f1_flat_signed = _flat([Rq.centered(c) for c in f1])
        if max(abs(c) for c in f1_flat_signed) > par.B_f1:
            return False
        z_b_signed = [Rqh.centered(c) for c in z_b]
        if max(abs(c) for c in _flat(z_b_signed)) > 6 * par.phi_b * par.B_b:
            return False

        # -- line 5: norm bounds on z_tilde, r_tilde, e_tilde.
        # Use the triangle bound for e_u = z''_u + r''_u: width
        # sigma_0 + sigma_0_prime.  This is intentionally looser than the
        # compact max-width expression and admits honest signatures with
        # comfortable slack.
        t = par.tail_t
        B_z = par.sigma_0 * t * math.sqrt(par.T * par.d * par.l)
        B_r = par.sigma_0_prime * t * math.sqrt(
            par.T * par.d * par.l_prime)
        B_e = (par.sigma_0 + par.sigma_0_prime) * t * math.sqrt(
            par.T * par.d * par.k)

        if math.sqrt(Rq.vec_l2_norm_sq(z_tilde)) > B_z:
            return False
        if math.sqrt(Rq.vec_l2_norm_sq(r_tilde)) > B_r:
            return False
        if math.sqrt(Rq.vec_l2_norm_sq(e_tilde)) > B_e:
            return False

        # -- lines 7-9: reconstruct  f_{j,0}  from  f_{j,1..beta-1}
        # f_full is in R_q for the selector polynomial computation.
        f_full = self._reconstruct_f(x, f1)

        # f0 infinity-norm check (paper Fig. 6 line 4).  Uses the tail
        # bound  B_f0 = 1 + tau_f0 * sqrt(beta-1) * sigma_a.
        f0_flat_signed = _flat(
            [Rq.centered(f_full[j][0]) for j in range(par.kappa)])
        if max(abs(c) for c in f0_flat_signed) > par.B_f0:
            return False

        # -- lines 10-12: compute  g_{j,i}  in R_qhat and check bounds
        # g values belong to the binary proof (R_qhat).
        g0, g1 = self._compute_g_qhat(x, f_full)
        if Rqh.vec_inf_norm(g0) > par.B_g0:
            return False
        if Rqh.vec_inf_norm(g1) > par.B_g1:
            return False

        # -- line 15: reconstruct  A_hat_bin  from binary proof
        g_all = self._interleave_g(g0, g1)
        f_all = self._flatten_f(f_full)

        # f values are small (bounded by norm check above) so
        # embedding from signed integers into R_qhat is exact.
        vec_for_G = Rqh.vec_concat(
            z_b,
            [Rqh.from_centered(Rq.centered(c)) for c in f_all],
            g_all,
        )
        Gv = Rqh.mat_vec(G, vec_for_G)

        # subtract  x * 2^{K_B} * B_bin^(1)
        x_hat = Rqh.from_centered(Rq.centered(x))
        scale_val = (1 << par.K_B)
        A_hat_bin = []
        for i in range(par.n_hat):
            t = Rqh.sub(Gv[i],
                        Rqh.scale(scale_val,
                                  Rqh.mul(x_hat, B_bin_hi[i])))
            A_hat_bin.append(t)

        # -- line 16: decompose A_hat_bin
        A_hat_hi = []
        for comp in A_hat_bin:
            h, lo = Rqh.centered_decompose(comp, par.K_A)
            # check |low| bound  (line 37 of Sign_bin)
            lo_signed = Rqh.centered(lo)
            bound = (1 << (par.K_A - 1)) - par.w * (1 << (par.K_B - 1))
            if max(abs(c) for c in lo_signed) > bound:
                return False
            A_hat_hi.append(h)

        # -- line 19: B := H_com(PK, mu)
        B = self._hash_com(pk_hash, mu)
        B_bar = self._augment_I(B, par.k)

        # -- line 20: aggregated keys
        agg_keys = self.kagg(pk_table, pk_hash)

        A = self._expand_A(pp)
        A_bar = self._augment_I(A, par.k)

        # lhs = A_bar z_tilde + B_bar r_tilde + e_tilde
        full_z = z_tilde + [Rq.zero()] * par.k
        full_r = r_tilde + [Rq.zero()] * par.k
        lhs_A = Rq.mat_vec(A_bar, full_z)
        lhs_B = Rq.mat_vec(B_bar, full_r)
        lhs = Rq.vec_add(lhs_A, lhs_B)
        lhs = Rq.vec_add(lhs, e_tilde)

        # -- line 21: reconstruct  w_hat_0  [new paper, Fig. 6]
        #   w_hat_0 = sum_i (prod_j f_{j,i_j}) t_tilde_i
        #             - (A_bar z_tilde + B_bar r_tilde + e_tilde)
        #             - sum_{j>=1} x^j 2^{K_{w,j}} w_tilde_j^(1)
        x_pow = [Rq.one()]
        for _ in range(kappa):
            x_pow.append(Rq.mul(x_pow[-1], x))

        pk_sum = Rq.vec_zero(par.k)
        for i in range(par.N):
            digits = self._base_digits(i, par.beta, par.kappa)
            selector = Rq.one()
            for j in range(kappa):
                selector = Rq.mul(selector, f_full[j][digits[j]])
            pk_sum = Rq.vec_add(pk_sum,
                                Rq.vec_scale(selector, agg_keys[i]))

        w_hat_0 = Rq.vec_sub(pk_sum, lhs)
        for j in range(1, kappa):
            scale_kw = 1 << par.K_w
            term = Rq.vec_scale_int(
                scale_kw,
                Rq.vec_scale(x_pow[j], w_tilde_hi_rest[j - 1]))
            w_hat_0 = Rq.vec_sub(w_hat_0, term)

        # -- line 22: decompose  w_hat_0 = 2^{K_{w,0}} w_hat_0^(1) + w_hat_0^(0)
        w_hat_0_hi = []
        for comp in w_hat_0:
            h, _ = Rq.centered_decompose(comp, par.K_w)
            w_hat_0_hi.append(h)

        # -- line 23-24: Fiat-Shamir check with reconstructed w_hat_0^(1)
        w_hi_for_hash = [w_hat_0_hi] + list(w_tilde_hi_rest)
        x_check = self._hash_fs(
            mu, A_hat_hi, B_bin_hi, w_hi_for_hash, pk_hash)
        if x_check != x:
            return False

        return True

    # ==================================================================
    #  Convenience  --  full two-round signing ceremony
    # ==================================================================

    def sign(self, pp, sks, ell, mu, pk_table, signing_seed):
        """
        Run the complete interactive signing protocol.

        sks[u]  = secret key of signer u  (row u of column ell).
        Returns the aggregated signature  sigma_tilde, or raises on
        exhausting max_attempts.
        """
        par = self.par
        T = par.T
        assert len(sks) == T
        pk_hash = self._pk_hash(pk_table)

        for attempt in range(par.max_attempts):
            rho = derive_subseed(signing_seed, b"rho", attempt)

            # ---- round 1 ----
            states = []
            all_coms = []
            for u in range(T):
                st, com = self.sign1(pp, sks[u], u, ell, mu,
                                     pk_table, rho, attempt, pk_hash)
                states.append(st)
                all_coms.append(com)

            # ---- round 2 ----
            sigmas = []
            restart = False
            for u in range(T):
                sig_u = self.sign2(states[u], all_coms, pk_table)
                if sig_u is None:
                    restart = True
                    break
                sigmas.append(sig_u)

            if restart:
                continue

            # ---- aggregate ----
            return self.sagg(sigmas)

        raise RuntimeError(
            f"signing failed after {par.max_attempts} attempts")

    # ==================================================================
    #  Sign_bin  (binary selection proof)           [Fig. 5, right]
    # ==================================================================

    def _sign_bin(self, pp, ell, mu, w_tilde_hi, pk_table,
                  rho, attempt, pk_hash=None):
        """
        Sign_bin(pp, ell, mu, (w_tilde_j^(1)), PK; rho).
        G is expanded from rho_pp; no separate rho_G is needed.

        Returns the proof transcript  pi  (dict), or None on Restart.
        """
        par = self.par
        Rq, Rqh = self.Rq, self.Rqh
        kappa, beta, d = par.kappa, par.beta, par.d
        if pk_hash is None:
            pk_hash = self._pk_hash(pk_table)

        # -- line 2: N = beta^kappa
        N = par.N

        # -- line 3: G := Expand(rho_pp, "G"; ...)
        G = self._expand_G(pp)

        # -- lines 4-6: decompose ell, form one-hot  b
        ell_digits = self._base_digits(ell, beta, kappa)
        b_vecs = []          # b[j][i]  for j in [0,kappa), i in [0,beta)
        for j in range(kappa):
            row = [Rqh.zero() for _ in range(beta)]
            row[ell_digits[j]] = Rqh.one()
            b_vecs.append(row)

        # -- line 7: b_1  (non-trivial entries  i in [1, beta-1])
        b1_flat = []
        for j in range(kappa):
            for i in range(1, beta):
                b1_flat.append(b_vecs[j][i])

        # -- line 8: expand  a_{j,i}  for  i in [1, beta-1]
        a_coeffs = self._expand_a_coeffs(rho, attempt)
        # a_coeffs[j][i] for j in [0,kappa), i in [0,beta)
        # (a_{j,0} already set to  -sum a_{j,i>0})

        # -- lines 12-14: form  a, c, d  vectors  (over R_qhat)
        a_flat = []    # all a_{j,i}  j in [0,kappa), i in [0,beta)
        c_flat = []    # a_{j,i} * (1 - 2*b_{j,i})
        d_flat = []    # -(a_{j,i})^2
        for j in range(kappa):
            for i in range(beta):
                a_ji = Rqh.from_centered(a_coeffs[j][i])
                a_flat.append(a_ji)
                if b_vecs[j][i] == Rqh.one():
                    # 1 - 2*1 = -1
                    c_flat.append(Rqh.neg(a_ji))
                else:
                    c_flat.append(a_ji)
                d_flat.append(Rqh.neg(Rqh.mul(a_ji, a_ji)))

        # -- line 15: r_b <- S_1^{n_hat+k_hat}  (ternary)
        xof_rb = make_xof(rho, b"rb", attempt)
        r_b = [Rqh.from_centered(
                   xof_sample_bounded(xof_rb, 1, d))
               for _ in range(par.n_hat + par.k_hat)]

        # -- line 15 cont: r_a <- D_{sigma_b}^{n_hat+k_hat}
        xof_ra = make_xof(rho, b"ra", attempt)
        r_a = [Rqh.from_centered(
                   xof_sample_gaussian(xof_ra, self.cdt_b, par.lam, d))
               for _ in range(par.n_hat + par.k_hat)]

        # -- line 16: B_bin := G (r_b, b, c)^T
        b_flat = []
        for j in range(kappa):
            for i in range(beta):
                b_flat.append(b_vecs[j][i])
        vec_commit = Rqh.vec_concat(r_b, b_flat, c_flat)
        B_bin = Rqh.mat_vec(G, vec_commit)

        # -- line 17: decompose  B_bin
        B_bin_hi, B_bin_lo = [], []
        for comp in B_bin:
            h, lo = Rqh.centered_decompose(comp, par.K_B)
            B_bin_hi.append(h)
            B_bin_lo.append(lo)

        # -- line 18: A_bin := G (r_a, a, d)^T
        vec_a = Rqh.vec_concat(r_a, a_flat, d_flat)
        A_bin = Rqh.mat_vec(G, vec_a)

        # -- line 19: decompose A_bin
        A_bin_hi = []
        for comp in A_bin:
            h, _ = Rqh.centered_decompose(comp, par.K_A)
            A_bin_hi.append(h)

        # -- line 20: Fiat-Shamir challenge
        x = self._hash_fs(mu, A_bin_hi, B_bin_hi, w_tilde_hi,
                                  pk_hash)

        # -- line 21: z_b := r_a + x * r_b
        x_hat = Rqh.from_centered(self.Rq.centered(x))
        z_b = Rqh.vec_add(r_a, Rqh.vec_scale(x_hat, r_b))

        # -- line 22: RejOp on z_b
        z_b_signed = [Rqh.centered(c) for c in z_b]
        v_signed = [Rqh.centered(c)
                    for c in Rqh.vec_scale(x_hat, r_b)]
        xof_rej_b = make_xof(rho, b"rej_b", attempt)
        if not rej_op(xof_rej_b,
                      _flat(z_b_signed), _flat(v_signed),
                      par.phi_b, par.B_b):
            return None

        # -- lines 23-28: f_{j,i} = x * delta_{ell_j, i} + a_{j,i}
        #
        # f values are small (bounded by 6*phi_a*B_a + 1) so they
        # are valid in both R_q and R_qhat without reduction.  We
        # store them as *signed integer lists* and convert to the
        # appropriate ring only when needed.
        f_full = []   # f_full[j][i] = signed coefficient list
        x_signed = self.Rq.centered(x)
        for j in range(kappa):
            row = []
            for i in range(beta):
                f_ji = list(a_coeffs[j][i])   # copy of signed list
                if i == ell_digits[j]:
                    f_ji = [f_ji[c] + x_signed[c] for c in range(d)]
                row.append(f_ji)
            f_full.append(row)

        # f1_flat for Rej: signed lists for i >= 1
        f1_signed = []
        for j in range(kappa):
            for i in range(1, beta):
                f1_signed.append(f_full[j][i])

        # -- line 29: Rej on f1
        b1_signed = []   # x * b_1  (shift for f_1)
        for j in range(kappa):
            for i in range(1, beta):
                if i == ell_digits[j]:
                    b1_signed.append(list(x_signed))
                else:
                    b1_signed.append([0] * d)
        xof_rej_a = make_xof(rho, b"rej_a", attempt)
        if not rej(xof_rej_a,
                   _flat(f1_signed),
                   _flat(b1_signed),
                   par.phi_a, par.B_a):
            return None

        # Infinity-norm restart checks on f_1 and f_0 (paper Fig. 5
        # Sign_bin line 12).  In an honest execution both bounds are
        # satisfied except with probability eps_each = eps_tot/4, but we
        # enforce them explicitly so the final signature matches Vf.
        f0_signed = [f_full[j][0] for j in range(kappa)]
        if max(abs(c) for c in _flat(f1_signed)) > par.B_f1:
            return None
        if max(abs(c) for c in _flat(f0_signed)) > par.B_f0:
            return None

        # -- lines 30-33: g_{j,i} = f_{j,i} * (x - f_{j,i})
        #
        # g values are part of the binary proof and belong in R_qhat.
        # We compute the product in Z via _poly_mul_signed, then
        # reduce into R_qhat.
        g0_list, g1_list = [], []   # R_qhat elements
        for j in range(kappa):
            for i in range(beta):
                f_ji = f_full[j][i]
                diff = [x_signed[c] - f_ji[c] for c in range(d)]
                g_ji = _poly_mul_signed(f_ji, diff, d)
                g_ji_reduced = Rqh.from_centered(g_ji)
                if i == 0:
                    g0_list.append(g_ji_reduced)
                else:
                    g1_list.append(g_ji_reduced)

        if Rqh.vec_inf_norm(g0_list) > par.B_g0:
            return None
        if Rqh.vec_inf_norm(g1_list) > par.B_g1:
            return None

        # -- lines 35-37: A_hat_bin check
        # f values → R_qhat (small, so from_centered is exact)
        f_all = []
        for j in range(kappa):
            for i in range(beta):
                f_all.append(Rqh.from_centered(f_full[j][i]))
        # g values are already in R_qhat
        g_all = []
        for j in range(kappa):
            g_all.append(g0_list[j])
            for i in range(1, beta):
                idx = j * (beta - 1) + (i - 1)
                g_all.append(g1_list[idx])

        vec_for_G = Rqh.vec_concat(z_b, f_all, g_all)
        Gv = Rqh.mat_vec(G, vec_for_G)

        A_hat_bin = []
        scale_kb = 1 << par.K_B
        for i in range(par.n_hat):
            t = Rqh.sub(Gv[i],
                        Rqh.scale(scale_kb,
                                  Rqh.mul(x_hat, B_bin_hi[i])))
            A_hat_bin.append(t)

        for comp in A_hat_bin:
            _, lo = Rqh.centered_decompose(comp, par.K_A)
            lo_s = Rqh.centered(lo)
            bound = (1 << (par.K_A - 1)) - par.w * (1 << (par.K_B - 1))
            if max(abs(c) for c in lo_s) > bound:
                return None

        # -- line 38: return pi
        # New paper: w_tilde_hi[0] is NOT included in the signature;
        # it is reconstructed by the verifier.  For kappa=1 this means
        # no w_tilde_hi terms are transmitted.
        #
        # f1 is stored as R_q elements (used by verifier for selector
        # polynomial in R_q; also embedded into R_qhat for the binary
        # proof check -- the values are small enough for both rings).
        f1_rq = [Rq.from_centered(s) for s in f1_signed]
        return dict(
            B_bin_hi=B_bin_hi,
            w_tilde_hi=w_tilde_hi[1:],
            x=x,
            f1=f1_rq,
            z_b=z_b,
        )

    # ==================================================================
    #  Internal helpers
    # ==================================================================

    def _expand_A(self, pp):
        """Expand public matrix  A in R_q^{k x l}  from seed."""
        par, Rq = self.par, self.Rq
        A = []
        for i in range(par.k):
            row = []
            for j in range(par.l):
                xof = make_xof(pp, b"A", i, j)
                row.append(xof_sample_uniform(xof, par.q, par.d))
            A.append(row)
        return A

    def _expand_G(self, pp):
        """Expand binary-proof matrix  G in R_{q_hat}^{n_hat x G_cols}."""
        par, Rqh = self.par, self.Rqh
        G = []
        for i in range(par.n_hat):
            row = []
            for j in range(par.G_cols):
                xof = make_xof(pp, b"G", i, j)
                row.append(xof_sample_uniform(xof, par.q_hat, par.d))
            G.append(row)
        return G

    def _augment_I(self, M, k):
        """
        [M | I_k]  -- append a  k x k  identity block on the right.

        M is  k  rows;  each row gets  k  extra ring elements.
        """
        Rq = self.Rq
        out = []
        for i, row in enumerate(M):
            new_row = list(row)
            for j in range(k):
                new_row.append(Rq.one() if i == j else Rq.zero())
            out.append(new_row)
        return out

    def _hash_agg(self, pk_hash, u):
        """
        H_agg(H(PK), u) -> challenge-set element alpha_u in R_q.
        """
        h = make_xof(pk_hash, b"agg", u)
        raw = xof_sample_challenge(h, self.par.w, self.par.d)
        return self.Rq.from_centered(raw)

    def _hash_com(self, pk_hash, mu):
        """
        H_com(H(PK), mu) -> B in R_q^{k x l'}.
        """
        par, Rq = self.par, self.Rq
        seed = derive_subseed(pk_hash, b"com", mu)
        B = []
        for i in range(par.k):
            row = []
            for j in range(par.l_prime):
                xof = make_xof(seed, b"B", i, j)
                row.append(xof_sample_uniform(xof, par.q, par.d))
            B.append(row)
        return B

    def _hash_fs(self, mu, A_hi, B_hi, w_hi, pk_hash):
        """
        Fiat-Shamir hash -> challenge  x in C.

        H(mu, A_bin^(1), B_bin^(1), w_hat_0^(1), w_tilde_1^(1), ..., H(PK))

        Coefficients are serialised as 8-byte signed LE to accommodate
        both R_q and R_qhat moduli up to ~2^63.
        """
        Rq = self.Rq
        h = make_xof(b"", b"FS")
        h.update(mu if isinstance(mu, bytes) else mu.encode())
        for comp in A_hi:
            for c in comp:
                h.update(int(c).to_bytes(8, "little", signed=True))
        for comp in B_hi:
            for c in comp:
                h.update(int(c).to_bytes(8, "little", signed=True))
        for j_vec in w_hi:
            for comp in j_vec:
                for c in comp:
                    h.update(int(c).to_bytes(8, "little", signed=True))
        h.update(pk_hash)
        raw = xof_sample_challenge(h, self.par.w, self.par.d)
        return self.Rq.from_centered(raw)

    def _pk_hash(self, pk_table):
        """
        256-bit digest of the canonical PK-table serialization.
        """
        return derive_subseed(self._pk_table_bytes(pk_table), b"pk")

    def _pk_table_bytes(self, pk_table):
        """Deterministic serialisation of PK for hashing.

        Uses 8-byte unsigned LE per coefficient to support moduli up
        to ~2^64.
        """
        parts = []
        for col in pk_table:
            for pk in col:
                for poly in pk:
                    for c in poly:
                        parts.append(int(c).to_bytes(8, "little"))
        return b"".join(parts)

    def _expand_a_coeffs(self, rho, attempt):
        """
        Expand Gaussian mask coefficients  a_{j,i}.

        Returns a_coeffs[j][i] as *signed* Python lists  (not mod q).
        For each j, a_{j,0} = -sum_{i>=1} a_{j,i}.
        """
        par = self.par
        kappa, beta, d = par.kappa, par.beta, par.d
        a_coeffs = []
        for j in range(kappa):
            row = [None] * beta
            for i in range(1, beta):
                xof = make_xof(rho, b"a", j, i, attempt)
                row[i] = xof_sample_gaussian(xof, self.cdt_a,
                                             par.lam, d)
            # a_{j,0} := -sum_{i>=1} a_{j,i}
            acc = [0] * d
            for i in range(1, beta):
                for c in range(d):
                    acc[c] += row[i][c]
            row[0] = [-acc[c] for c in range(d)]
            a_coeffs.append(row)
        return a_coeffs

    def _compute_p_coeffs(self, a_coeffs, ell):
        """
        Selector polynomial coefficients  p_{i,j}.

        For kappa = 1:  p_i(x) = prod_{j=0}^{kappa-1} (a_{j,i_j} + X b_{j,i_j})
        where the product is just one factor, so p_{i,0} is the constant
        term and the degree-kappa coefficient tells us if  i == ell.

        Returns  p_coeffs[i][j]  for  i in [0,N), j in [0,kappa)  as
        ring elements in R_q (unsigned canonical).
        """
        par = self.par
        Rq, d = self.Rq, self.par.d
        kappa, beta = par.kappa, par.beta
        N = par.N

        ell_digits = self._base_digits(ell, beta, kappa)

        p_coeffs = []
        for i in range(N):
            i_digits = self._base_digits(i, beta, kappa)
            # For kappa = 1 the "product" is a single linear factor:
            #   p_i(X) = a_{0, i_0}  +  X * delta_{i_0, ell_0}
            # Coefficients in x^0 .. x^{kappa-1}:
            coeffs_j = []
            if kappa == 1:
                coeffs_j.append(
                    Rq.from_centered(a_coeffs[0][i_digits[0]]))
            else:
                # General kappa: expand the product symbolically.
                # p_i(X) = prod_{j=0}^{kappa-1} (a_{j,i_j} + X b_{j,i_j})
                # where b_{j,i_j} = delta_{i_j, ell_j}
                # This is a degree-kappa polynomial in X; we only need
                # coefficients 0 .. kappa-1.
                poly_x = [[0] * d for _ in range(kappa + 1)]
                poly_x[0] = list(a_coeffs[0][i_digits[0]])
                b_val = 1 if i_digits[0] == ell_digits[0] else 0
                if b_val:
                    poly_x[1] = [b_val] * 1 + [0] * (d - 1)
                    # Actually b is a scalar (0 or 1), constant poly
                    poly_x[1] = [b_val] + [0] * (d - 1)
                for j_idx in range(1, kappa):
                    a_new = a_coeffs[j_idx][i_digits[j_idx]]
                    b_new = (1 if i_digits[j_idx] == ell_digits[j_idx]
                             else 0)
                    new_poly = [[0] * d for _ in range(kappa + 1)]
                    for deg in range(kappa + 1):
                        if any(c != 0 for c in poly_x[deg]):
                            # multiply by a_new  ->  same degree
                            prod_a = _poly_mul_signed(
                                poly_x[deg], a_new, d)
                            for c in range(d):
                                new_poly[deg][c] += prod_a[c]
                            # multiply by X*b_new  ->  degree+1
                            if b_new and deg + 1 <= kappa:
                                for c in range(d):
                                    new_poly[deg + 1][c] += \
                                        poly_x[deg][c] * b_new
                    poly_x = new_poly
                for j_idx in range(kappa):
                    coeffs_j.append(Rq.from_centered(poly_x[j_idx]))
            p_coeffs.append(coeffs_j)
        return p_coeffs

    def _reconstruct_f(self, x, f1_flat):
        """
        Reconstruct full  f_{j,i}  from  f_1  (the  i>=1  entries)
        and the challenge x.

        f_{j,0} = x - sum_{i>=1} f_{j,i}     [Vf line 8]

        Returns  f_full[j][i]  as ring elements (unsigned canonical).
        """
        par, Rq = self.par, self.Rq
        kappa, beta = par.kappa, par.beta

        idx = 0
        f_full = []
        for j in range(kappa):
            row = [None] * beta
            acc = Rq.zero()
            for i in range(1, beta):
                row[i] = f1_flat[idx]
                acc = Rq.add(acc, row[i])
                idx += 1
            row[0] = Rq.sub(x, acc)      # f_{j,0} = x - sum f_{j,i>0}
            f_full.append(row)
        return f_full

    def _compute_g(self, x, f_full):
        """g_{j,i} = f_{j,i} * (x - f_{j,i})  in R_q.  Returns (g0, g1)."""
        par, Rq = self.par, self.Rq
        g0, g1 = [], []
        for j in range(par.kappa):
            for i in range(par.beta):
                diff = Rq.sub(x, f_full[j][i])
                g = Rq.mul(f_full[j][i], diff)
                if i == 0:
                    g0.append(g)
                else:
                    g1.append(g)
        return g0, g1

    def _compute_g_qhat(self, x, f_full):
        """
        g_{j,i} = f_{j,i} * (x - f_{j,i})  in R_qhat.

        f_full contains R_q elements but the values are small (norm-
        checked by the verifier).  We convert to signed integers,
        multiply in Z, and reduce mod q_hat.  Returns (g0, g1) as
        R_qhat elements.
        """
        par = self.par
        Rq, Rqh = self.Rq, self.Rqh
        d = par.d
        x_signed = Rq.centered(x)
        g0, g1 = [], []
        for j in range(par.kappa):
            for i in range(par.beta):
                f_signed = Rq.centered(f_full[j][i])
                diff = [x_signed[c] - f_signed[c] for c in range(d)]
                g_signed = _poly_mul_signed(f_signed, diff, d)
                g_qhat = Rqh.from_centered(g_signed)
                if i == 0:
                    g0.append(g_qhat)
                else:
                    g1.append(g_qhat)
        return g0, g1

    def _interleave_g(self, g0, g1):
        """Flatten g0, g1 into the order expected by G multiplication."""
        par = self.par
        out = []
        idx1 = 0
        for j in range(par.kappa):
            out.append(g0[j])
            for _ in range(1, par.beta):
                out.append(g1[idx1])
                idx1 += 1
        return out

    def _flatten_f(self, f_full):
        """Flatten f_full[j][i] into a single list."""
        out = []
        for j_row in f_full:
            out.extend(j_row)
        return out

    @staticmethod
    def _base_digits(n, base, num_digits):
        """Decompose n in base 'base', least-significant first."""
        digits = []
        for _ in range(num_digits):
            digits.append(n % base)
            n //= base
        return digits


# ---- module-level helpers ------------------------------------------------

def _poly_mul_signed(a, b, d):
    """
    Negacyclic convolution of two *signed* coefficient lists.

    Returns a signed list (not reduced mod q).  Used for computing
    selector polynomials and quadratic terms over Z.
    """
    c = [0] * d
    for i in range(d):
        ai = a[i]
        if ai == 0:
            continue
        for j in range(d):
            k = i + j
            if k < d:
                c[k] += ai * b[j]
            else:
                c[k - d] -= ai * b[j]
    return c


# --------------------------------------------------------------------------
if __name__ == "__main__":
    from params import TEST_PARAMS

    print("Initialising LoTRS (building CDT tables) ...")
    scheme = LoTRS(TEST_PARAMS)

    pp = scheme.setup(b"\x00" * 32)

    # key generation for a small ring
    N, T = TEST_PARAMS.N, TEST_PARAMS.T
    pk_table = []           # pk_table[col][row] = pk
    all_sks = []            # all_sks[col][row] = sk
    for col in range(N):
        col_pks, col_sks = [], []
        for row in range(T):
            seed = bytes([col, row]) + b"\x00" * 30
            sk, pk = scheme.keygen(pp, seed)
            col_sks.append(sk)
            col_pks.append(pk)
        pk_table.append(col_pks)
        all_sks.append(col_sks)

    # sign with column 1
    ell = 1
    sks = all_sks[ell]
    mu = b"hello LoTRS"

    print("Signing ...")
    sig = scheme.sign(pp, sks, ell, mu, pk_table,
                      b"\xAA" * 32)
    print("Verifying ...")
    ok = scheme.verify(pp, mu, sig, pk_table)
    print(f"Verification: {'PASS' if ok else 'FAIL'}")
    assert ok, "honest signature must verify"

    print("lotrs.py: all self-tests passed")
