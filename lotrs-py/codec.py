"""
codec.py -- Compact serialization for LoTRS objects.

Fixed-width packing for uniform / bit-dropped components,
Golomb-Rice coding for Gaussian-distributed components.

All decoders enforce canonical form: out-of-range coefficients,
nonzero padding bits, and trailing bytes are rejected with
ValueError so that data from untrusted channels cannot produce
surprising internal state.
"""

import math

from ring import Ring
from params import LoTRSParams


# =====================================================================
#  Bitstream helpers
# =====================================================================

class BitWriter:
    """Accumulates bits (LSB-first within each byte) and flushes to bytes."""

    __slots__ = ("buf",)

    def __init__(self):
        self.buf = []

    def write_bits(self, value, n):
        """Write *n* bits of *value*, least-significant first."""
        for i in range(n):
            self.buf.append((value >> i) & 1)

    def write_unary(self, value):
        """Unary code: *value* one-bits followed by a zero-bit."""
        for _ in range(value):
            self.buf.append(1)
        self.buf.append(0)

    def pad_to_byte(self):
        while len(self.buf) % 8:
            self.buf.append(0)

    def to_bytes(self):
        self.pad_to_byte()
        out = bytearray(len(self.buf) // 8)
        for i, b in enumerate(self.buf):
            if b:
                out[i >> 3] |= 1 << (i & 7)
        return bytes(out)

    def bit_length(self):
        return len(self.buf)


class BitReader:
    """Reads bits (LSB-first) from a bytes object."""

    __slots__ = ("data", "pos")

    def __init__(self, data):
        self.data = data
        self.pos = 0                       # bit position

    def read_bits(self, n):
        value = 0
        for i in range(n):
            byte_idx = self.pos >> 3
            bit_idx = self.pos & 7
            if byte_idx >= len(self.data):
                raise ValueError("truncated data")
            value |= ((self.data[byte_idx] >> bit_idx) & 1) << i
            self.pos += 1
        return value

    def read_unary(self, max_val):
        """Read unary-coded value; reject if > *max_val* (DoS guard)."""
        count = 0
        while True:
            b = self.read_bits(1)
            if b == 0:
                return count
            count += 1
            if count > max_val:
                raise ValueError(
                    f"unary run {count} exceeds limit {max_val}")

    def check_padding(self):
        """Advance to next byte boundary; reject nonzero padding."""
        while self.pos & 7:
            if self.read_bits(1) != 0:
                raise ValueError("nonzero padding bits")

    def check_exhausted(self):
        """Reject if any bytes remain unread."""
        self.check_padding()
        if (self.pos >> 3) != len(self.data):
            raise ValueError(
                f"{len(self.data) - (self.pos >> 3)} trailing bytes")

    @property
    def consumed_bytes(self):
        return (self.pos + 7) >> 3


# =====================================================================
#  Fixed-width polynomial packing
# =====================================================================

def _pack_fixed(coeffs, dx, offset=0):
    """
    Pack *d* integer coefficients at *dx* bits each.

    Before packing, each coefficient is shifted by +offset so that
    negative values become non-negative (caller sets offset = bound).
    """
    w = BitWriter()
    for c in coeffs:
        v = c + offset
        if v < 0 or v >= (1 << dx):
            raise ValueError(
                f"coefficient {c} (shifted {v}) out of {dx}-bit range")
        w.write_bits(v, dx)
    w.pad_to_byte()
    return w.to_bytes()


def _unpack_fixed(data, d, dx, offset=0):
    """
    Unpack *d* coefficients at *dx* bits each.

    Returns (list[int], bytes_consumed).  Rejects nonzero padding.
    """
    r = BitReader(data)
    coeffs = []
    for _ in range(d):
        v = r.read_bits(dx)
        coeffs.append(v - offset)
    r.check_padding()
    return coeffs, r.consumed_bytes


# =====================================================================
#  Golomb-Rice polynomial packing
# =====================================================================

def optimal_rice_k(sigma):
    """Optimal Rice parameter for discrete Gaussian with std dev *sigma*."""
    if sigma < 1:
        return 0
    return max(0, int(math.floor(math.log2(1.1774 * sigma))))


def _pack_rice(coeffs, rice_k, bound):
    """
    Golomb-Rice encode a signed-integer polynomial.

    For each coefficient c:
      abs_c = |c|
      low   = abs_c  &  ((1 << rice_k) - 1)
      high  = abs_c  >>  rice_k
      bits  = low (rice_k bits) + unary(high) + sign (if c != 0)

    Byte-aligned after all d coefficients.
    """
    w = BitWriter()
    for c in coeffs:
        abs_c = abs(c)
        if abs_c > bound:
            raise ValueError(
                f"|{c}| = {abs_c} exceeds Rice bound {bound}")
        low = abs_c & ((1 << rice_k) - 1)
        high = abs_c >> rice_k
        w.write_bits(low, rice_k)
        w.write_unary(high)
        if abs_c != 0:
            w.write_bits(1 if c < 0 else 0, 1)
    w.pad_to_byte()
    return w.to_bytes()


def _unpack_rice(data, d, rice_k, bound):
    """
    Golomb-Rice decode *d* signed-integer coefficients.

    Returns (list[int], bytes_consumed).  Rejects out-of-range and
    nonzero padding.
    """
    max_high = (bound >> rice_k) + 1       # DoS guard
    r = BitReader(data)
    coeffs = []
    for _ in range(d):
        low = r.read_bits(rice_k)
        high = r.read_unary(max_high)
        abs_c = (high << rice_k) | low
        if abs_c > bound:
            raise ValueError(
                f"Rice-decoded |coeff| {abs_c} exceeds bound {bound}")
        if abs_c == 0:
            coeffs.append(0)
        else:
            sign = r.read_bits(1)
            coeffs.append(-abs_c if sign else abs_c)
    r.check_padding()
    return coeffs, r.consumed_bytes


# =====================================================================
#  Challenge encoding
# =====================================================================

def _pack_challenge(poly, w, d):
    """
    Encode a challenge polynomial (w nonzero +/-1 coefficients).

    Format: w positions (sorted ascending, ceil(log2(d)) bits each)
            + w sign bits (0 = +1, 1 = -1), byte-aligned.
    """
    pos_bits = max(1, (d - 1).bit_length())
    positions = sorted(i for i, c in enumerate(poly) if c != 0)
    if len(positions) != w:
        raise ValueError(
            f"challenge has {len(positions)} nonzero coeffs, expected {w}")

    wr = BitWriter()
    for p in positions:
        wr.write_bits(p, pos_bits)
    for p in positions:
        wr.write_bits(1 if poly[p] < 0 else 0, 1)
    wr.pad_to_byte()
    return wr.to_bytes()


def _unpack_challenge(data, w, d):
    """
    Decode a challenge polynomial.

    Rejects duplicate / out-of-order positions, non-binary values.
    Returns (list[int], bytes_consumed).
    """
    pos_bits = max(1, (d - 1).bit_length())
    r = BitReader(data)

    positions = []
    for _ in range(w):
        p = r.read_bits(pos_bits)
        if p >= d:
            raise ValueError(f"challenge position {p} >= d={d}")
        if positions and p <= positions[-1]:
            raise ValueError("challenge positions not strictly increasing")
        positions.append(p)

    signs = []
    for _ in range(w):
        signs.append(r.read_bits(1))

    r.check_padding()

    poly = [0] * d
    for p, s in zip(positions, signs):
        poly[p] = -1 if s else 1

    return poly, r.consumed_bytes


# =====================================================================
#  LoTRS codec
# =====================================================================

class LoTRSCodec:
    """
    Serialize / deserialize all LoTRS objects.

    Encoding parameters are derived from the scheme's LoTRSParams.
    """

    def __init__(self, par: LoTRSParams):
        self.par = par
        self.Rq = Ring(par.q, par.d)
        self.Rqh = Ring(par.q_hat, par.d)

        # fixed-width bit widths
        self.dx_pk = (par.q - 1).bit_length()          # pk coefficients
        self.dx_bbin = (par.q_hat - 1).bit_length() - par.K_B
        self.dx_whi = (par.q - 1).bit_length() - par.K_w

        # Rice parameters and bounds for Gaussian components
        self.rice_f1 = optimal_rice_k(par.sigma_a)
        self.bound_f1 = int(math.ceil(6 * par.sigma_a))

        self.rice_zb = optimal_rice_k(par.sigma_b)
        self.bound_zb = int(math.ceil(6 * par.sigma_b))

        # e_tilde has effective width sigma_0 + sigma_0_prime
        # (triangle inequality on z''_u + r''_u summed over u);
        # see the matching note in lotrs.py verify().
        sigma_z = par.sigma_0 * math.sqrt(par.T)
        sigma_r = par.sigma_0_prime * math.sqrt(par.T)
        sigma_e = (par.sigma_0 + par.sigma_0_prime) * math.sqrt(par.T)

        t = par.tail_t
        self.rice_zt = optimal_rice_k(sigma_z)
        self.bound_zt = int(math.ceil(
            t * par.sigma_0 * math.sqrt(par.T * par.d * par.l)))

        self.rice_rt = optimal_rice_k(sigma_r)
        self.bound_rt = int(math.ceil(
            t * par.sigma_0_prime
            * math.sqrt(par.T * par.d * par.l_prime)))

        self.rice_et = optimal_rice_k(sigma_e)
        self.bound_et = int(math.ceil(
            t * (par.sigma_0 + par.sigma_0_prime)
            * math.sqrt(par.T * par.d * par.k)))

    # ---- public parameters -----------------------------------------------

    def pp_encode(self, pp):
        """Encode public parameters (32-byte seed)."""
        if len(pp) != 32:
            raise ValueError("pp must be 32 bytes")
        return bytes(pp)

    def pp_decode(self, data):
        """Decode public parameters."""
        if len(data) != 32:
            raise ValueError(f"pp must be 32 bytes, got {len(data)}")
        return bytes(data)

    # ---- secret key (seed form) ------------------------------------------

    def sk_encode(self, seed):
        if len(seed) != 32:
            raise ValueError("sk seed must be 32 bytes")
        return bytes(seed)

    def sk_decode(self, data):
        if len(data) != 32:
            raise ValueError(f"sk must be 32 bytes, got {len(data)}")
        return bytes(data)

    # ---- public key ------------------------------------------------------

    def pk_encode(self, pk):
        """Encode public key (k ring elements in R_q)."""
        parts = []
        for poly in pk:
            parts.append(_pack_fixed(poly, self.dx_pk))
        return b"".join(parts)

    def pk_decode(self, data):
        """Decode and validate public key."""
        par = self.par
        poly_bytes = (self.dx_pk * par.d + 7) // 8
        expected = par.k * poly_bytes
        if len(data) != expected:
            raise ValueError(
                f"pk: expected {expected} bytes, got {len(data)}")
        pk = []
        for i in range(par.k):
            chunk = data[i * poly_bytes:(i + 1) * poly_bytes]
            coeffs, _ = _unpack_fixed(chunk, par.d, self.dx_pk)
            # validate range
            for c in coeffs:
                if c < 0 or c >= par.q:
                    raise ValueError(f"pk coeff {c} out of [0, {par.q})")
            pk.append(coeffs)
        return pk

    # ---- signature -------------------------------------------------------

    def sig_encode(self, sigma):
        """
        Encode an aggregated LoTRS signature.

        Layout (all byte-aligned per component):
          B_bin_hi    : n_hat polys, fixed-width (dx_bbin)
          w_tilde_hi  : kappa * k polys, fixed-width (dx_whi)
          x           : challenge encoding
          f1          : kappa*(beta-1) polys, Rice
          z_b         : (n_hat+k_hat) polys, Rice
          z_tilde     : l polys, Rice
          r_tilde     : l' polys, Rice
          e_tilde     : k polys, Rice
        """
        par = self.par
        pi = sigma["pi"]
        parts = []

        # B_bin^(1) -- fixed width over R_{q_hat}
        for poly in pi["B_bin_hi"]:
            centered = self.Rqh.centered(poly)
            offset = (1 << (self.dx_bbin - 1))
            parts.append(_pack_fixed(centered, self.dx_bbin, offset))

        # w_tilde^(1) -- fixed width over R_q
        for j_vec in pi["w_tilde_hi"]:
            for poly in j_vec:
                centered = self.Rq.centered(poly)
                offset = (1 << (self.dx_whi - 1))
                parts.append(
                    _pack_fixed(centered, self.dx_whi, offset))

        # x -- challenge
        parts.append(_pack_challenge(
            self.Rq.centered(pi["x"]), par.w, par.d))

        # f1 -- Rice
        for poly in pi["f1"]:
            centered = self.Rq.centered(poly)
            parts.append(
                _pack_rice(centered, self.rice_f1, self.bound_f1))

        # z_b -- Rice
        for poly in pi["z_b"]:
            centered = self.Rqh.centered(poly)
            parts.append(
                _pack_rice(centered, self.rice_zb, self.bound_zb))

        # z_tilde -- Rice
        for poly in sigma["z_tilde"]:
            centered = self.Rq.centered(poly)
            parts.append(
                _pack_rice(centered, self.rice_zt, self.bound_zt))

        # r_tilde -- Rice
        for poly in sigma["r_tilde"]:
            centered = self.Rq.centered(poly)
            parts.append(
                _pack_rice(centered, self.rice_rt, self.bound_rt))

        # e_tilde -- Rice
        for poly in sigma["e_tilde"]:
            centered = self.Rq.centered(poly)
            parts.append(
                _pack_rice(centered, self.rice_et, self.bound_et))

        return b"".join(parts)

    def sig_decode(self, data):
        """
        Decode and validate an aggregated LoTRS signature.

        Raises ValueError on any malformation.
        """
        par = self.par
        Rq, Rqh = self.Rq, self.Rqh
        pos = 0

        def consume_fixed(ring, n_polys, dx, q):
            nonlocal pos
            poly_bytes = (dx * par.d + 7) // 8
            polys = []
            offset = (1 << (dx - 1))
            for _ in range(n_polys):
                chunk = data[pos:pos + poly_bytes]
                if len(chunk) < poly_bytes:
                    raise ValueError("truncated signature")
                coeffs, _ = _unpack_fixed(chunk, par.d, dx, offset)
                # reduce to unsigned canonical in [0, q)
                polys.append([c % q for c in coeffs])
                pos += poly_bytes
            return polys

        def consume_rice(ring, n_polys, rice_k, bound, q):
            nonlocal pos
            polys = []
            for _ in range(n_polys):
                coeffs, consumed = _unpack_rice(
                    data[pos:], par.d, rice_k, bound)
                polys.append(ring.from_centered(coeffs))
                pos += consumed
            return polys

        # B_bin^(1)
        B_bin_hi = consume_fixed(
            Rqh, par.n_hat, self.dx_bbin, par.q_hat)

        # w_tilde^(1) -- only j >= 1 entries (kappa-1 vectors)
        n_whi = (par.kappa - 1) * par.k
        w_hi_flat = consume_fixed(Rq, n_whi, self.dx_whi, par.q)
        w_tilde_hi = []
        idx = 0
        for _ in range(par.kappa - 1):
            w_tilde_hi.append(w_hi_flat[idx:idx + par.k])
            idx += par.k

        # x (challenge)
        x_centered, consumed = _unpack_challenge(
            data[pos:], par.w, par.d)
        x = Rq.from_centered(x_centered)
        pos += consumed

        # f1
        f1 = consume_rice(Rq, par.kappa * (par.beta - 1),
                          self.rice_f1, self.bound_f1, par.q)

        # z_b
        z_b = consume_rice(Rqh, par.n_hat + par.k_hat,
                           self.rice_zb, self.bound_zb, par.q_hat)

        # z_tilde
        z_tilde = consume_rice(Rq, par.l,
                               self.rice_zt, self.bound_zt, par.q)

        # r_tilde
        r_tilde = consume_rice(Rq, par.l_prime,
                               self.rice_rt, self.bound_rt, par.q)

        # e_tilde
        e_tilde = consume_rice(Rq, par.k,
                               self.rice_et, self.bound_et, par.q)

        if pos != len(data):
            raise ValueError(
                f"{len(data) - pos} trailing bytes in signature")

        return dict(
            pi=dict(
                B_bin_hi=B_bin_hi,
                w_tilde_hi=w_tilde_hi,
                x=x,
                f1=f1,
                z_b=z_b,
            ),
            z_tilde=z_tilde,
            r_tilde=r_tilde,
            e_tilde=e_tilde,
        )

    # ---- size reporting --------------------------------------------------

    def sizes(self):
        """Estimated byte sizes for each signature component."""
        par = self.par
        d = par.d

        def rice_est(n_polys, rice_k, sigma):
            # average bits per coefficient ≈ rice_k + 1 + 1/mean_high + sign
            # ≈ log2(4.13 * sigma) for optimal rice_k
            avg_bits = math.log2(4.13 * max(sigma, 1)) + 1
            return int(math.ceil(n_polys * d * avg_bits / 8))

        sigma_z = par.sigma_0 * math.sqrt(par.T)
        sigma_r = par.sigma_0_prime * math.sqrt(par.T)
        sigma_e = (par.sigma_0 + par.sigma_0_prime) * math.sqrt(par.T)

        pos_bits = max(1, (d - 1).bit_length())
        ch_bits = par.w * (pos_bits + 1)

        return {
            "B_bin_hi":    par.n_hat * (self.dx_bbin * d + 7) // 8,
            "w_tilde_hi":  (par.kappa - 1) * par.k
                           * (self.dx_whi * d + 7) // 8,
            "x":           (ch_bits + 7) // 8,
            "f1":          rice_est(par.kappa * (par.beta - 1),
                                   self.rice_f1, par.sigma_a),
            "z_b":         rice_est(par.n_hat + par.k_hat,
                                   self.rice_zb, par.sigma_b),
            "z_tilde":     rice_est(par.l, self.rice_zt, sigma_z),
            "r_tilde":     rice_est(par.l_prime, self.rice_rt, sigma_r),
            "e_tilde":     rice_est(par.k, self.rice_et, sigma_e),
        }

    def total_size_estimate(self):
        s = self.sizes()
        return sum(s.values())

    def print_sizes(self):
        s = self.sizes()
        total = sum(s.values())
        print(f"{'Component':14s}  {'Bytes':>7s}  {'%':>5s}  Encoding")
        print("-" * 52)
        enc = {
            "B_bin_hi":   f"fixed {self.dx_bbin}b",
            "w_tilde_hi": f"fixed {self.dx_whi}b",
            "x":          f"challenge w={self.par.w}",
            "f1":         f"Rice k={self.rice_f1}",
            "z_b":        f"Rice k={self.rice_zb}",
            "z_tilde":    f"Rice k={self.rice_zt}",
            "r_tilde":    f"Rice k={self.rice_rt}",
            "e_tilde":    f"Rice k={self.rice_et}",
        }
        for name in ["B_bin_hi", "w_tilde_hi", "x", "f1", "z_b",
                      "z_tilde", "r_tilde", "e_tilde"]:
            b = s[name]
            print(f"  {name:12s}  {b:7d}  {100*b/total:5.1f}  {enc[name]}")
        print("-" * 52)
        print(f"  {'TOTAL':12s}  {total:7d}")
        return total


# --------------------------------------------------------------------------
if __name__ == "__main__":
    from params import TEST_PARAMS

    codec = LoTRSCodec(TEST_PARAMS)

    print("Signature size estimate (TEST_PARAMS):")
    codec.print_sizes()

    # quick round-trip on Rice primitives
    import random
    rng = random.Random(42)
    sigma = 100.0
    rk = optimal_rice_k(sigma)
    bound = int(6 * sigma)
    coeffs = [int(rng.gauss(0, sigma)) for _ in range(32)]
    coeffs = [max(-bound, min(bound, c)) for c in coeffs]
    enc = _pack_rice(coeffs, rk, bound)
    dec, _ = _unpack_rice(enc, 32, rk, bound)
    assert dec == coeffs, "Rice round-trip"

    # challenge round-trip
    ch = [0] * 32
    ch[3] = 1
    ch[7] = -1
    ch[15] = 1
    ch[28] = -1
    enc_ch = _pack_challenge(ch, 4, 32)
    dec_ch, _ = _unpack_challenge(enc_ch, 4, 32)
    assert dec_ch == ch, "challenge round-trip"

    print("\ncodec.py: all self-tests passed")
