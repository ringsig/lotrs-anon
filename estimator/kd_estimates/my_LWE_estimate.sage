load("MLWE_security.py")

d = 256			# Ring dimension of Rq
q = 2^27.5		# modulus
k = 4			# secret key dimension over Rq
m = 4			# number of samples over Rq
B = 1			# error is sampled from {-B,..,B}




print("*********** Running MLWE estimates ************")

params = MLWEParameterSet(d, k, m, B, q, distr="uniform")
bkz_block = MLWE_summarize_attacks(params)[0]
print("Root Hermite Factor =", round(delta_BKZ(bkz_block),5))
print()
