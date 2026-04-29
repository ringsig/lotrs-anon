from estimator import *
from sage.all import *

from kd_estimates.MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from kd_estimates.model_BKZ import delta_BKZ
from ASIS_sec_estimate.ASIS_MSIS_security import MSIS_summarize_attacks as ASIS_summarize_attacks
from ASIS_sec_estimate.ASIS_MSIS_security import MSISParameterSet as ASISParameterSet
from ASIS_sec_estimate.ASIS_model_BKZ import delta_BKZ as ASIS_delta_BKZ

from lotrs_finder import calculate_PK, calculate_sig_size, number_reps, setBinASISBounds, setDualMSASISBounds

from lotrs_param_checks import *


"""Auxiliary functions"""
def prime_5_mod_8(bits):
    q = int(next_prime(int(2) ** int(bits - 1)))
    while q % 8 != 5:
        q = int(next_prime(q))
    return q

"""Main function"""
def main():
    RHF_max = 1.0045

    """Parameter variables"""
    #Common parameters
    T = 50 #Threshold size
    N = 100 #Ring size; n = beta^kappa
    kappa = 1
    beta = N #kappa = 1 --> beta = N

    t=1.2 #Gaussian tail bound param

    B=1 #All MLWE instances use B=1

    x_inf_norm = 1 #Infinity norm of secret and error in MLWE instances
    w = 31 #Hamming weight of challenge

    eta_s, eta_prime_s = 1, 1 #Smoothing param

    d = 128 #dim(R_q), i.e., R_q = Z_q[X]/(X^d+1)

    size_x = ceil(log(binomial(d,w),2)+w) #Number of bits for challenge

    #Binary proof dimensions + parameters
    nhat = 10 #sisrank
    khat = 8
    if False:
        approx_logqhat = 34
        qhat = prime_5_mod_8(approx_logqhat)
    qhat = 8589934237
    logq_hat = RR(log(qhat, 2))
    print("Binary proof MLWE modulus q_hat =", qhat, "logq_hat", logq_hat)
    print("nhat =", nhat, "khat =", khat)

    phi_a = 50
    phi_b = 4

    #Dropped bits
    K_b = 5
    K_w0 = 5
    K_a = round(log((nhat*d*w*pow(2, K_b)), 2))
    print("K_a" , K_a)

    #Set probability that bounds on f_0, f_1, g_0, g_1 will reject
    eps_total = RR(0.01)

    #Binary proof bounds
    bin_arr = setBinASISBounds(beta, kappa, d, w, nhat, khat, phi_a, phi_b, K_b, eps_total)
    B1 = bin_arr[0][0]
    m1 = bin_arr[0][1]

    B2 = bin_arr[1][0]
    m2 = bin_arr[1][1]

    B3 = bin_arr[3][0] #We will skip the third value, i.e. index 2
    m3 = bin_arr[3][1]

    B4 = bin_arr[4][0]
    m4 = bin_arr[4][1]

    B5 = bin_arr[5][0]
    m5 = bin_arr[5][1]

    #Binary std deviations
    B_b = sqrt(d* nhat * khat)*w
    B_a = sqrt(kappa*w)
    sigma_zb = phi_b*B_b
    sigma_f = phi_a*B_a

    #DualMS dimensions + parameters
    l = 5
    l_prime = l+1
    k = 12
    if False:
        approx_logq = 39
        q_dualms = prime_5_mod_8(approx_logq)
    q_dualms = 274877906837
    logq = RR(log(q_dualms, 2))
    print("\nDualMS modulus q =", q_dualms, "logq = ", logq)
    print("l =", l, "l_prime =", l_prime, "k =", k)

    phi = 11.75 * T

    b2KB = pow(2,13) #Convert bits to bytes

    #Check whether to compress \tilde{w}_0
    compress_w = True

    #DualMS bounds
    sigma_s =  ceil(RR(((2*d)/sqrt(2*pi))*pow(q_dualms, (k/(l+k)+2/(d*(l+k))))))
    sigma_s_prime = ceil(RR(((2*d)/sqrt(2*pi))*pow(q_dualms, (k/(l_prime+k)+2/(d*(l_prime+k))))))

    B_zero = RR(sqrt(d*(l+k))*(eta_s*pow(w, kappa+1)+6*sigma_s*(kappa-1)*pow(w, kappa-1)))
    B_zero_prime = RR(sqrt(d*(l_prime+k))*(eta_prime_s*pow(w, kappa+1)+6*sigma_s_prime*(kappa-1)*pow(w, kappa-1)))
    B_hat_0 = max(B_zero, B_zero_prime)

    B_z = RR(t*phi*B_zero*sqrt(d*l))
    B_tilde_z = RR(sqrt(T)*B_z)
    #print("B_tilde_z =", log(pow(B_tilde_z, 2), 2))

    B_r = RR(t*phi*B_zero_prime*sqrt(d*l_prime))
    B_tilde_r = RR(sqrt(T)*B_r)
    #print("B_tilde_r =", log(pow(B_tilde_r, 2), 2))

    B_e = RR(t*phi*B_hat_0*sqrt(d*k))
    B_tilde_e = RR(sqrt(T)*B_e)
    #print("B_tilde_e =", log(pow(B_tilde_e, 2), 2))

    B_det = RR(pow((2*w), kappa*(kappa+1)/2))
    print("log2(B_det) =", log(pow(B_det, 2), 2))

    B_Gamma = RR(pow((2*w), kappa*(kappa-1)/2))

    #MSIS bound for DualMS MSIS instance
    #Compression residual bound for kappa = 1 case only
    if(compress_w):
        B_eta_w = RR(pow(2, K_w0-1))
        print("log2(k*B_eta_w) =", log((k*pow(B_eta_w,2)), 2)) #Need to compare other terms with k*pow(B_eta_w,2)
        beta_sis = RR(sqrt(4*d*pow(B_det, 2)+ 4*d*pow((kappa+1),2)*(pow(B_Gamma,2)*(pow(B_tilde_z, 2)+pow(B_tilde_r, 2)+pow(B_tilde_e, 2)+k*pow(B_eta_w, 2)))))
        total_width = 1+l+l_prime+2*k
    else:
        #k*pow(B_eta_w, 2) term missing when no compression on \tilde{w}_j0 done
        beta_sis = RR(sqrt(4*d*pow(B_det, 2)+ 4*d*pow((kappa+1),2)*(pow(B_Gamma,2)*(pow(B_tilde_z, 2)+pow(B_tilde_r, 2)+pow(B_tilde_e, 2)))))
        total_width = 1+l+l_prime+k

    #DualMS MLWE std deviation
    sigma_b0 = phi*B_zero
    sigma_tilde_z = phi*B_zero*sqrt(T)
    sigma_tilde_r = phi*B_zero_prime*sqrt(T)
    sigma_tilde_e = phi*B_hat_0*sqrt(T)

    print("\n=== Binary Proof LWE Estimator ===")
    bin_lwe_params = LWE.Parameters(n=(d*khat), q=qhat, Xs=ND.Uniform(-1,1), Xe=ND.DiscreteGaussian(sigma_zb))
    print(bin_lwe_params)
    results_bin_lwe = LWE.estimate.rough(bin_lwe_params)
    print(results_bin_lwe)

    print("\n\n=== Binary Proof ASIS Estimator ===")
    attack_variant = 0
    for attack_variant in [0,1,2]:
        # function below assumes B1>=B2>=B3>=B4>=B5
        print("***** Attack Variant", attack_variant, " *********")
        params = ASISParameterSet(d, m1+m2+m3+m4+m5, nhat, B1, B2, B3, B4, B5, m1, m2, m3, m4, m5, qhat, norm="linf")
        (m_pq, b_pq, c_pq) = ASIS_summarize_attacks(params, attack_variant=attack_variant)
        print("BKZ block size =", b_pq)
        print("RHF =", round(ASIS_delta_BKZ(b_pq),5))
        print("Cost=", round(c_pq,5))
        print()


    print("\n=== DualMS LWE Estimator ===")
    dualms_lwe_params = LWE.Parameters(n=(d*l), q=q_dualms, Xs=ND.Uniform(-1,1), Xe=ND.DiscreteGaussian(sigma_b0))
    print(dualms_lwe_params)
    results_dualms_params = LWE.estimate.rough(dualms_lwe_params)
    print(results_dualms_params)

    # print("\n\n=== DualMS MSIS Estimator ===")
    # try:
    #     print("Running MSIS_summarize_attacks w.r.t. l2 norm...")
    #     params = MSISParameterSet(d, total_width, k, beta_sis, q_dualms, norm="l2")
    #     bkz_block = MSIS_summarize_attacks(params)[0]
    #     RHF = round(delta_BKZ(bkz_block),5)
    #     print("Root Hermite Factor =", RHF)
    # except Exception as e:
    #     print("MSIS_summarize_attacks failed:", e)

    print("\n\n=== DualMS ASIS Estimator ===")
    arrDualMS = setDualMSASISBounds(T, kappa, d, w, k, l, l_prime, eta_s, eta_prime_s, t, logq, phi, K_w0, compress_w)
    B1, m1 = arrDualMS[0]
    B2, m2 = arrDualMS[1]
    B3, m3 = arrDualMS[2]
    B4, m4 = arrDualMS[3]
    B5, m5 = arrDualMS[4]

    for attack_variant in [0,1,2]:
        # function below assumes B1>=B2>=B3>=B4>=B5
        print("***** Attack Variant", attack_variant, " *********")
        params = ASISParameterSet(d, m1+m2+m3+m4+m5, k, B1, B2, B3, B4, B5, m1, m2, m3, m4, m5, q_dualms, norm="linf")
        (m_pq, b_pq, c_pq) = ASIS_summarize_attacks(params, attack_variant=attack_variant)
        print("BKZ block size =", b_pq)
        print("RHF =", round(ASIS_delta_BKZ(b_pq),5))
        print("Cost=", round(c_pq,5))
        print()
    
    sig_size = calculate_sig_size(kappa, beta, nhat, khat, k, l, l_prime, d, logq_hat, logq, K_b, K_w0, size_x, sigma_f, sigma_zb, sigma_tilde_z, sigma_tilde_r, sigma_tilde_e, compress_w)
    print("\nSignature size:", RR(sig_size/b2KB), "KB\n")

    size_single_pk, size_PK = calculate_PK(T, N, k, d, logq)
    print("Single public key size:", RR(size_single_pk/b2KB), "KB\n")
    print("Ring PK size:", RR(size_PK/b2KB), "KB\n")

    number_lotrs_reps = number_reps(T, phi_a, phi_b, phi, nhat, d, w, K_a, K_b, K_w0, eps_total, compress_w)
    print("\nNumber of repetitions for rejection sampling:", number_lotrs_reps)

    print("\n\n=== Perform condition checks ===")
    print("q_hat:", check_q_prime_5_mod_8(qhat), "\n")
    print("q:", check_q_prime_5_mod_8(q_dualms), "\n")
    checkChallengeDiff(qhat, q_dualms, x_inf_norm)
    print("Regularity condition for sigma_s:", check_sigma(sigma_s, d, q_dualms, k, l))
    print("Regularity condition for sigma_s_prime:", check_sigma(sigma_s_prime, d, q_dualms, k, l_prime), "\n")
    checkRangeProofCondition(d, kappa, phi_a, w, qhat)

    

"""------------------------------------------------------------------------------------------------------------------"""
if __name__ == "__main__":
    main()