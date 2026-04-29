from sage.all import *

"""
1. q ≡ 5 (mod 8)  ⇒  r = 2 and q ≡ 1 (mod 2r)
2. Check Δ is invertible in R_q, R_qhat
3. Check ||Δ·b||_∞ ≤ 2B
4. Check q > d(α||Δ||_∞ + 2B)^2 (range proof condition from MatRiCT+)
   - usually take α = 1, so q > 2(||Δ||_∞ + 2B)^2
5. Check regularity condition holds
"""

def check_q_prime_5_mod_8(q):
    rem = mod(q,8) #rem should be a multiple of 5
    if(rem == 5):
        print(q, " equiv 5 mod 8")
        return True
    else:
        return False

def check_sigma(sigma_s, d, q, k, l):
    exponent = k/(l+k)+2/(d*(l+k))
    sigma_check = RR((2*d/sqrt(2*pi))*pow(q, exponent))
    return (sigma_s > sigma_check)

def checkChallengeDiff(qhat, q, x_inf_norm):
    print("Challenge diffs invertible mod qhat:", RR(2*x_inf_norm) < sqrt(RR(qhat)/ 2))
    print("Challenge diffs invertible mod q", RR(2*x_inf_norm) < sqrt(RR(q)/ 2), "\n")

def checkRangeProofCondition(d, kappa, phi_a, w, qhat):
    B_f1 = RR(6 * phi_a * sqrt(kappa * w))
    qhat_min = (d*pow(2+2*B_f1, 2))
    print("log(qhat_min)",log(qhat_min, 2).n())
    print("log(qhat)",log(qhat, 2).n())
    print("Range proof condition:", qhat > qhat_min)