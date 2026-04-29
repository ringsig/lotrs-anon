"""
sample.py -- Deterministic XOF-based samplers for LoTRS.

Every sampler consumes bytes from a caller-supplied SHAKE XOF, so all
outputs are fully reproducible given the XOF seed.  This enables
known-answer-test (KAT) generation.
"""

import math
from dataclasses import dataclass
from Crypto.Hash import SHAKE256

try:
    import mpmath
    _HAS_MPMATH = True
except ImportError:
    _HAS_MPMATH = False


# ---- CDT construction ----------------------------------------------------

DEFAULT_GAUSSIAN_TAILCUT = 14
MAX_CDT_ENTRIES = 1_000_000
FACCT_U_BITS = 48
FACCT_PROB_BITS = 64
FACCT_EXP_POLY_DEGREE = 20
FACCT_EXP_POLY_Q64 = [
    18446744073709551616,
    -18446744073709551616,
    9223372036854775808,
    -3074457345618258432,
    768614336404564608,
    -153722867280912928,
    25620477880152156,
    -3660068268593165,
    457508533574146,
    -50834281508238,
    5083428150824,
    -462129831893,
    38510819324,
    -2962370717,
    211597908,
    -14106527,
    881658,
    -51862,
    2881,
    -152,
    8,
]
SHAKE256_RATE = 136
SHAKE128_RATE = 168
DEFAULT_XOF_BATCH_BLOCKS = 4


@dataclass(frozen=True)
class GaussianSampler:
    """
    Prepared Gaussian sampler description.

    kind = "cdt"
        data is a CDT list as used by xof_sample_gaussian().

    kind = "facct"
        data is a dict of fixed-point constants for the large-sigma
        FACCT-style sampler.  This is only selected for the masking
        widths sigma_0 / sigma_0' when a direct table would be
        impractically large.
    """
    kind: str
    data: object
    sigma: float
    lam: int


class BufferedXOF:
    """
    Buffered wrapper around a SHAKE XOF.

    PyCryptodome already squeezes in rate-sized blocks internally. This
    wrapper batches small reads in multiples of the permutation rate so
    sampler hot loops avoid repeated tiny Python/C calls while keeping
    the exact same output byte stream.
    """
    __slots__ = ("_shake", "_rate", "_batch_blocks", "_buf", "_pos")

    def __init__(self, shake, rate, batch_blocks=DEFAULT_XOF_BATCH_BLOCKS):
        self._shake = shake
        self._rate = rate
        self._batch_blocks = max(1, batch_blocks)
        self._buf = b""
        self._pos = 0

    def _refill(self, needed):
        unread = len(self._buf) - self._pos
        if unread:
            self._buf = self._buf[self._pos:]
        else:
            self._buf = b""
        self._pos = 0
        fetch = max(self._batch_blocks * self._rate,
                    ((needed + self._rate - 1) // self._rate) * self._rate)
        self._buf += self._shake.read(fetch)

    def read(self, n):
        if n < 0:
            raise ValueError("read length must be non-negative")
        if n == 0:
            return b""
        unread = len(self._buf) - self._pos
        if unread < n:
            self._refill(n)
        out = self._buf[self._pos:self._pos + n]
        self._pos += n
        return out

    def update(self, data):
        """
        Absorb more input before any squeezing has happened.
        """
        if self._buf or self._pos:
            raise RuntimeError("cannot update XOF after read()")
        self._shake.update(data)
        return self

def build_cdt(sigma, lam=128, tailcut=14):
    """
    Cumulative-distribution table for  |D_sigma|.

    cdt[k] = floor( Pr[|X| <= k] * 2^lam )   for k = 0 .. tailcut*sigma.

    Last entry is clamped to 2^lam so binary search always terminates.
    Requires *mpmath* for the high-precision CDF computation.
    """
    if not _HAS_MPMATH:
        raise RuntimeError("mpmath is required for CDT construction")

    prec_digits = max(2 * lam + 64, 256) // 3 + 10
    bound = int(math.ceil(tailcut * sigma))

    with mpmath.workdps(prec_digits):
        two_s2 = 2 * mpmath.mpf(sigma) ** 2
        scale = mpmath.power(2, lam)

        # un-normalised probabilities  rho(k) = exp(-k^2 / 2 sigma^2)
        rho = [mpmath.exp(-mpmath.mpf(k * k) / two_s2)
               for k in range(bound + 1)]
        total = rho[0] + 2 * mpmath.fsum(rho[1:])

        # running CDF
        cdf = []
        running = rho[0]
        for k in range(bound + 1):
            if k > 0:
                running += 2 * rho[k]
            cdf.append(int(running * scale / total))

        cdf[-1] = 1 << lam          # termination guarantee
    return cdf


def needs_large_gaussian_sampler(sigma, tailcut=DEFAULT_GAUSSIAN_TAILCUT,
                                 max_entries=MAX_CDT_ENTRIES):
    """
    Predicate: would a direct CDT for this sigma exceed *max_entries*?

    This is **informational** — it documents the cutoff the scheme-level
    code uses to pick between CDT and FACCT.  Unlike the old
    `prepare_gaussian_sampler`, no sampler builder branches on this
    silently; callers decide explicitly (see `build_cdt_sampler` and
    `build_facct_sampler` below).
    """
    return int(math.ceil(tailcut * sigma)) + 1 > max_entries


def _validate_sigma(sigma):
    sigma = float(sigma)
    if not math.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("sigma must be a finite positive float")
    return sigma


def build_cdt_sampler(sigma, lam=128, tailcut=DEFAULT_GAUSSIAN_TAILCUT):
    """
    Build an exact-CDT Gaussian sampler for this sigma.

    Caller responsibility: this requires `ceil(14*sigma) + 1` table
    entries in memory.  Use `needs_large_gaussian_sampler(sigma)` to
    check feasibility, or just call `build_facct_sampler` for the
    masking widths `sigma_0` / `sigma_0_prime` at BENCH / PRODUCTION.

    This is the only path that requires `mpmath`.
    """
    sigma = _validate_sigma(sigma)
    return GaussianSampler("cdt",
                           build_cdt(sigma, lam=lam, tailcut=tailcut),
                           sigma, lam)


def build_facct_sampler(sigma, lam=128, tailcut=DEFAULT_GAUSSIAN_TAILCUT):
    """
    Build a FACCT-style large-sigma Gaussian sampler for this sigma.

    No CDT.  Runtime uses only integer arithmetic; prep builds a small
    fixed-point state dict (the polynomial coefficients themselves are
    frozen at `FACCT_EXP_POLY_Q64`).  See
    `../lotrs-facct-sampler.md` for the spec.
    """
    sigma = _validate_sigma(sigma)
    return GaussianSampler("facct",
                           _prepare_large_gaussian_sampler(
                               sigma, lam=lam, tailcut=tailcut),
                           sigma, lam)


# ---- XOF helpers ---------------------------------------------------------

def make_xof(seed, *tags):
    """
    Create a SHAKE256 instance with domain-separated seed.

    Tags may be  bytes | str | int.  Ints are encoded as 4-byte LE.
    """
    h = SHAKE256.new(seed)
    for tag in tags:
        if isinstance(tag, int):
            h.update(tag.to_bytes(4, "little"))
        elif isinstance(tag, str):
            h.update(tag.encode("ascii"))
        elif isinstance(tag, bytes):
            h.update(tag)
        else:
            raise TypeError(f"unsupported tag type {type(tag)}")
    return BufferedXOF(h, SHAKE256_RATE)


def derive_subseed(seed, *tags, outlen=32):
    """
    Derive a labeled sub-secret from a root seed via XOF expansion.

    This is a counter-/label-mode convenience helper:
    `seed || tags -> outlen bytes`.
    It is suitable for places that already want a single derived block
    and want explicit domain separation without sharing a mutable XOF.
    """
    return make_xof(seed, *tags).read(outlen)


# ---- polynomial samplers -------------------------------------------------

def xof_sample_uniform(xof, q, d):
    """Uniform polynomial in  [0, q)^d  via rejection."""
    bits_q = (q - 1).bit_length()
    byte_q = (bits_q + 7) // 8
    mask_q = (1 << bits_q) - 1

    coeffs = []
    while len(coeffs) < d:
        x = int.from_bytes(xof.read(byte_q), "little") & mask_q
        if x < q:
            coeffs.append(x)
    return coeffs


def xof_sample_gaussian(xof, cdt, lam, d):
    """
    Discrete Gaussian  D_sigma  via CDT look-up.

    For each coefficient: read *lam* bits of magnitude randomness,
    binary-search *cdt* for the magnitude, then read a separate sign
    byte.  The CDT is scaled to 2^lam, so the full lam-bit value is
    needed for an unbiased comparison.
    """
    lam_bytes = lam // 8
    coeffs = []
    for _ in range(d):
        u = int.from_bytes(xof.read(lam_bytes), "little")
        # separate sign byte: only the LSB is used
        sign = xof.read(1)[0] & 1

        # binary search: smallest k with  cdt[k] > u
        lo, hi = 0, len(cdt) - 1
        while lo < hi:
            mid = (lo + hi) >> 1
            if cdt[mid] <= u:
                lo = mid + 1
            else:
                hi = mid
        k = lo
        coeffs.append(-k if (k > 0 and sign) else k)
    return coeffs


def _sample_u01_from_xof(xof):
    """
    Sample u in (0, 1) from 53 mantissa bits of XOF output.

    Matches the same little-endian / 53-bit mantissa convention used by
    the rejection-sampling code below so the resulting sampler remains
    fully deterministic from the XOF stream.
    """
    u_raw = int.from_bytes(xof.read(8), "little")
    u = (u_raw & ((1 << 53) - 1)) / (1 << 53)
    if u == 0.0:
        u = 2.0 ** -53
    return u


def _xof_randbelow(xof, bound):
    """
    Uniform integer in [0, bound) from XOF bytes.
    """
    if bound <= 0:
        raise ValueError("bound must be positive")
    bits = max(1, (bound - 1).bit_length())
    nbytes = (bits + 7) // 8
    mask = (1 << bits) - 1
    while True:
        x = int.from_bytes(xof.read(nbytes), "little") & mask
        if x < bound:
            return x


def _prepare_large_gaussian_sampler(sigma, lam=128,
                                    tailcut=DEFAULT_GAUSSIAN_TAILCUT,
                                    u_bits=FACCT_U_BITS,
                                    prob_bits=FACCT_PROB_BITS,
                                    poly_degree=FACCT_EXP_POLY_DEGREE):
    """
    Prepare fixed-point constants for the large-sigma FACCT-style
    integer sampler.

    Runtime path:
    - propose a signed integer uniformly in [-ceil(tailcut*sigma), +ceil(tailcut*sigma)]
    - compute u = x^2 / (2 sigma^2)
    - decompose u = k ln 2 + r with r in [0, ln 2)
    - accept with probability 2^{-k} * exp(-r)

    The Bernoulli-exp remainder exp(-r) is approximated by a fixed-point
    Taylor polynomial on [0, ln 2).  All runtime operations are integer.
    """
    sigma_sq = sigma * sigma
    tail = int(math.ceil(tailcut * sigma))
    prob_scale = 1 << prob_bits
    ln2_q = int(round(math.log(2.0) * (1 << u_bits)))
    if prob_bits != FACCT_PROB_BITS:
        raise ValueError("FACCT exp polynomial is fixed for PROB_BITS = 64")
    if poly_degree != FACCT_EXP_POLY_DEGREE:
        raise ValueError("FACCT exp polynomial is fixed for degree 20")
    exp_poly = FACCT_EXP_POLY_Q64
    return dict(
        tail=tail,
        u_bits=u_bits,
        prob_bits=prob_bits,
        prob_scale=prob_scale,
        sigma_sq_int=int(round(sigma_sq)),
        ln2_q=ln2_q,
        exp_poly=exp_poly,
    )


def _facct_exp_neg_threshold(r_q, u_bits, prob_scale, exp_poly):
    """
    Approximate exp(-r) on r in [0, ln 2) using a fixed-point Taylor
    polynomial evaluated by Horner's rule.
    """
    acc = exp_poly[-1]
    for coeff in reversed(exp_poly[:-1]):
        acc = coeff + ((acc * r_q) >> u_bits)
    if acc < 0:
        return 0
    if acc >= prob_scale:
        return prob_scale - 1
    return acc


def _xof_bernoulli_power_of_two(xof, k):
    """
    Return True with probability 2^{-k} by requiring k XOF bits to be 0.
    """
    while k >= 64:
        if int.from_bytes(xof.read(8), "little") != 0:
            return False
        k -= 64
    if k == 0:
        return True
    nbytes = (k + 7) // 8
    mask = (1 << k) - 1
    return (int.from_bytes(xof.read(nbytes), "little") & mask) == 0


def xof_sample_gaussian_facct(xof, prepared, d):
    """
    Deterministic large-sigma FACCT-style sampler using integer arithmetic.

    This is the large-sigma path used for sigma_0 / sigma_0'.  It keeps
    the small exact CDT path untouched and uses only integer operations
    during runtime sampling.

    The target distribution is the truncated centered discrete Gaussian
    on |x| <= ceil(tailcut * sigma), with FACCT-style Bernoulli-exp
    decomposition for the acceptance probability.
    """
    tail = prepared["tail"]
    u_bits = prepared["u_bits"]
    prob_scale = prepared["prob_scale"]
    sigma_sq_int = prepared["sigma_sq_int"]
    ln2_q = prepared["ln2_q"]
    exp_poly = prepared["exp_poly"]
    coeffs = []
    while len(coeffs) < d:
        x = _xof_randbelow(xof, 2 * tail + 1) - tail
        u_q = ((x * x) << u_bits) // (2 * sigma_sq_int)
        k = u_q // ln2_q
        r_q = u_q - k * ln2_q
        if not _xof_bernoulli_power_of_two(xof, k):
            continue
        threshold = _facct_exp_neg_threshold(
            r_q, u_bits, prob_scale, exp_poly)
        if int.from_bytes(xof.read(8), "little") >= threshold:
            continue
        coeffs.append(x)
    return coeffs


def xof_sample_short(xof, eta, d):
    """Uniform polynomial with coefficients in  [-eta, eta]."""
    width = 2 * eta + 1
    bits = width.bit_length()
    nbytes = (bits + 7) // 8
    mask = (1 << bits) - 1

    coeffs = []
    while len(coeffs) < d:
        x = int.from_bytes(xof.read(nbytes), "little") & mask
        if x < width:
            coeffs.append(x - eta)
    return coeffs


def xof_sample_challenge(xof, w, d):
    """
    Challenge  x in C = { x in R : ||x||_inf = 1, ||x||_1 = w }.

    Picks *w* distinct positions, assigns each +/-1.

    Positions are sampled by rejection over the smallest power of 2
    that contains d, so there is no modulo bias even when d is not a
    power of 2.  For power-of-2 d the reject rate is zero.
    """
    pos_bits = max(1, (d - 1).bit_length())
    pos_bytes = (pos_bits + 7) // 8
    pos_mask = (1 << pos_bits) - 1

    positions = []
    seen = set()
    while len(positions) < w:
        raw = int.from_bytes(xof.read(pos_bytes), "little") & pos_mask
        if raw >= d:            # reject to avoid modulo bias
            continue
        if raw not in seen:
            seen.add(raw)
            positions.append(raw)

    sign_bytes = xof.read((w + 7) // 8)
    sign_bits = int.from_bytes(sign_bytes, "little")

    coeffs = [0] * d
    for i, pos in enumerate(positions):
        coeffs[pos] = -1 if ((sign_bits >> i) & 1) else 1
    return coeffs


def xof_sample_ternary(xof, d):
    """Uniform ternary polynomial (coefficients in {-1, 0, 1})."""
    coeffs = []
    while len(coeffs) < d:
        b = xof.read(1)[0]
        for shift in range(0, 8, 2):
            if len(coeffs) >= d:
                break
            val = (b >> shift) & 3
            if val < 3:                       # reject 3
                coeffs.append(val - 1)        # 0->-1, 1->0, 2->1
    return coeffs


def xof_sample_bounded(xof, bound, d):
    """Uniform polynomial with coefficients in  [-bound, bound]."""
    return xof_sample_short(xof, bound, d)


# ---- rejection sampling (Rej / RejOp, Fig. 1) ---------------------------

def _coeff_inner(a, b):
    """Inner product of two coefficient lists (plain integers)."""
    return sum(ai * bi for ai, bi in zip(a, b))


def _coeff_norm_sq(a):
    return sum(ai * ai for ai in a)


def _flat(vec_of_polys):
    """Flatten a vector of ring elements into a single coefficient list."""
    out = []
    for p in vec_of_polys:
        out.extend(p)
    return out


def rej(xof, z_flat, v_flat, phi, K):
    """
    Rejection sampling  Rej(z, v, phi, K)  -- Fig. 1, left column.

    Returns True  on accept (output 0 in the paper = continue),
            False on reject (output 1 = abort / restart).

    *z_flat* and *v_flat* are flattened signed-integer coefficient lists.
    """
    sigma = phi * K
    sigma_sq = sigma * sigma
    mu_phi = math.exp(12.0 / phi + 1.0 / (2.0 * phi * phi))

    inner_zv = _coeff_inner(z_flat, v_flat)
    norm_v_sq = _coeff_norm_sq(v_flat)

    log_p = (-2.0 * inner_zv + norm_v_sq) / (2.0 * sigma_sq)
    # p = (1 / mu_phi) * exp(log_p)
    # accept if  u < p   where u in [0, 1)
    log_threshold = log_p - math.log(mu_phi)

    # sample u in (0, 1)  from XOF
    u = _sample_u01_from_xof(xof)

    if math.log(u) >= log_threshold:
        return False                          # reject

    # infinity-norm check
    max_abs = max(abs(c) for c in z_flat) if z_flat else 0
    if max_abs > 6 * sigma:
        return False                          # reject

    return True                               # accept


def rej_op(xof, z_flat, v_flat, phi, K):
    """
    Optimised rejection sampling  RejOp(z, c, phi, K)  -- Fig. 1, right.

    Exploits  <z, v> >= 0.  Returns True on accept, False on reject.
    """
    inner_zv = _coeff_inner(z_flat, v_flat)
    if inner_zv < 0:
        return False                          # line 1: <z,v> < 0

    sigma = phi * K
    sigma_sq = sigma * sigma
    mu_phi = math.exp(1.0 / (2.0 * phi * phi))

    norm_v_sq = _coeff_norm_sq(v_flat)
    log_p = (-2.0 * inner_zv + norm_v_sq) / (2.0 * sigma_sq)
    log_threshold = log_p - math.log(mu_phi)

    u = _sample_u01_from_xof(xof)

    if math.log(u) >= log_threshold:
        return False

    max_abs = max(abs(c) for c in z_flat) if z_flat else 0
    if max_abs > 6 * sigma:
        return False

    return True


# --------------------------------------------------------------------------
if __name__ == "__main__":
    # determinism check
    xof1 = make_xof(b"\x00" * 32, b"test")
    xof2 = make_xof(b"\x00" * 32, b"test")
    assert xof_sample_uniform(xof1, 7681, 8) == \
           xof_sample_uniform(xof2, 7681, 8), "uniform determinism"

    # challenge has correct weight / norm
    xof3 = make_xof(b"\x01" * 32, b"chal")
    ch = xof_sample_challenge(xof3, 10, 64)
    assert sum(abs(c) for c in ch) == 10, "challenge L1"
    assert max(abs(c) for c in ch) == 1, "challenge Linf"

    # Gaussian has bounded coefficients
    cdt = build_cdt(10.0, lam=128)
    xof4 = make_xof(b"\x02" * 32, b"gauss")
    g = xof_sample_gaussian(xof4, cdt, 128, 256)
    assert max(abs(c) for c in g) <= 14 * 10, "Gaussian tail"

    # short polynomial bounds
    xof5 = make_xof(b"\x03" * 32, b"short")
    s = xof_sample_short(xof5, 3, 64)
    assert all(-3 <= c <= 3 for c in s), "short bounds"

    print("sample.py: all self-tests passed")
