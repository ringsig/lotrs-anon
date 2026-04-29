load("MSIS_security.py")


d = 32				# Ring dimension of Rq
w = 90				# width of the MSIS matrix (over Rq); i.e. number of columns
h = 53				# height of the MSIS matrix (over Rq); i.e. number of rows
B = 2^35			# *infinity* norm bound
B2 = 2^38			# *l2* norm bound
q = int(2^37.5)		# modulus q. q/2 needs to be larger than B. Typically q slightly larger than 2*B gives the best results for optimizing parameters
q2 = int(2^40.5)	# modulus for l2-norm-based check



print("*********** Running MSIS Estimates w.r.t. infinity norm ************")
params = MSISParameterSet(d, w, h, B, q, norm="linf")
bkz_block = MSIS_summarize_attacks(params)[0]
RHF = round(delta_BKZ(bkz_block),5)
print("Root Hermite Factor =", RHF)
print()


print("*********** Running MSIS Estimates w.r.t. l2 norm ************")
params = MSISParameterSet(d, w, h, B2, q2, norm="l2")
bkz_block = MSIS_summarize_attacks(params)[0]
RHF = round(delta_BKZ(bkz_block),5)
print("Root Hermite Factor =", RHF)
print()
