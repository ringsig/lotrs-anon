from kd_estimates.MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from kd_estimates.MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from math import sqrt


class UniformDilithiumParameterSet(object):
    def __init__(self, n, k, l, gamma, q, eta, pkdrop=0):
        self.n = n		# dim(Rq)
        self.k = k		# SIS rank
        self.l = l		# LWE rank
        self.gamma = gamma
        self.q = q
        self.eta = eta
        self.B = max(2*gamma, 2**(pkdrop+1))
        self.pkdrop = pkdrop


class GaussianDilithiumParameterSet(object):
    def __init__(self, n, k, l, sigma, q, eta, pkdrop=0):
        self.n = n
        self.k = k
        self.l = l
        self.sigma = sigma
        self.q = q
        self.eta = eta
        self.pkdrop = pkdrop
        self.B = 2*equation5(self)


def equation5(dps):
    B2 = ((1.05 * dps.sigma * sqrt((dps.k + dps.l)*dps.n))**2 
          +(2**(dps.pkdrop-1) * sqrt(60*dps.n*dps.k))**2)
    return sqrt(B2)

n = 256
q = 8380417
gamma = (q-1)/16

UnifWeakDilithium           = UniformDilithiumParameterSet(n, 3, 2, gamma, q, 7, pkdrop=14)
UnifMediumDilithium         = UniformDilithiumParameterSet(n, 4, 3, gamma, q, 6, pkdrop=14)
UnifRecommendedDilithium    = UniformDilithiumParameterSet(n, 5, 4, gamma, q, 5, pkdrop=14)
UnifVeryHighDilithium       = UniformDilithiumParameterSet(n, 6, 5, gamma, q, 3, pkdrop=14)


all_params_unif = [("Uniform Dilithium Weak", UnifWeakDilithium),
                   ("Uniform Dilithium Medium", UnifMediumDilithium),
                   ("Uniform Dilithium Recommended", UnifRecommendedDilithium),
                   ("Uniform Dilithium Very High", UnifVeryHighDilithium)]

all_params = all_params_unif

def Dilithium_to_MSIS(dps):
    if type(dps)==UniformDilithiumParameterSet:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.B, dps.q, norm="linf")
    if type(dps)==GaussianDilithiumParameterSet:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.B, dps.q, norm="l2")
    else:
        raise ValueError("Unrecognized Dilithium Parameter Type")


def Dilithium_to_MLWE(dps):
    return MLWEParameterSet(dps.n, dps.l, dps.k, dps.eta, dps.q, distr="uniform")

text_SIS = ["BKZ block-size $b$ to break SIS","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]
text_LWE = ["BKZ block-size $b$ to break LWE","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]



table_SIS = [4*[0] for i in range(4)]
table_LWE = [4*[0] for i in range(4)]
j = 0

for (scheme, param) in all_params_unif:
    print("\n"+scheme)
    print(param.__dict__)
    print("")
    v = MSIS_summarize_attacks(Dilithium_to_MSIS(param))
    for i in range(4):
        table_SIS[i][j] = v[i]
    v = MLWE_summarize_attacks(Dilithium_to_MLWE(param))
    for i in range(4):
        table_LWE[i][j] = v[i]
    j+=1

print("UNIFORM DILITHIUM TABLE")
print("========================")

print("\\hline")
for j in range(4):
    print(text_SIS[j]+" & "),
    for i in range(4):
        print(table_SIS[j][i]),
        if i<3:
            print(" & "),
    print("\\\\")
print("\\hline")
for j in range(4):
    print(text_LWE[j]+" & "),
    for i in range(4):
        print(table_LWE[j][i]),
        if i<3:
            print(" & "),
    print("\\\\")
print("\\hline")

print("========================")

table_SIS = [4*[0] for i in range(4)]
table_LWE = [4*[0] for i in range(4)]
j = 0



'''
************ OUTPUT ***************
Uniform Dilithium Weak
{'n': 256, 'k': 3, 'l': 2, 'gamma': 523776.0, 'q': 8380417, 'eta': 7, 'B': 1047552.0, 'pkdrop': 14}

Attack uses block-size 235 and 1536 dimensions, with 0 q-vectors
log2(epsilon) = -38.18, log2 nvector per run 48.77
shortest vector used has length l=16954523.19, q=8380417, `l<q'= 0
SIS & 1536 & 235 & 68 & 62 & 48
Primal attacks uses block-size 200 and 560 samples
Primal & 560 & 200 & 58 & 53 & 41
Dual attacks uses block-size 200 and 565 samples
shortest vector used has length l=1651835.59, q=8380417, `l<q'= 1
log2(epsilon) = -20.65, log2 nvector per run 41.50
Dual & 565 & 200 & 58 & 53 & 41

Uniform Dilithium Medium
{'n': 256, 'k': 4, 'l': 3, 'gamma': 523776.0, 'q': 8380417, 'eta': 6, 'B': 1047552.0, 'pkdrop': 14}

Attack uses block-size 355 and 2048 dimensions, with 0 q-vectors
log2(epsilon) = -48.06, log2 nvector per run 73.67
shortest vector used has length l=19353502.38, q=8380417, `l<q'= 0
SIS & 2048 & 355 & 103 & 94 & 73
Primal attacks uses block-size 345 and 765 samples
Primal & 765 & 345 & 100 & 91 & 71
Dual attacks uses block-size 345 and 755 samples
shortest vector used has length l=2502287.37, q=8380417, `l<q'= 1
log2(epsilon) = -35.54, log2 nvector per run 71.59
Dual & 755 & 345 & 100 & 91 & 71

Uniform Dilithium Recommended
{'n': 256, 'k': 5, 'l': 4, 'gamma': 523776.0, 'q': 8380417, 'eta': 5, 'B': 1047552.0, 'pkdrop': 14}

Attack uses block-size 475 and 2560 dimensions, with 0 q-vectors
log2(epsilon) = -92.38, log2 nvector per run 98.57
shortest vector used has length l=23121762.38, q=8380417, `l<q'= 0
SIS & 2560 & 475 & 138 & 125 & 98
Primal attacks uses block-size 485 and 1090 samples
Primal & 1090 & 485 & 141 & 128 & 100
Dual attacks uses block-size 485 and 1035 samples
shortest vector used has length l=3521975.56, q=8380417, `l<q'= 1
log2(epsilon) = -50.30, log2 nvector per run 100.65
Dual & 1035 & 485 & 141 & 128 & 100

Uniform Dilithium Very High
{'n': 256, 'k': 6, 'l': 5, 'gamma': 523776.0, 'q': 8380417, 'eta': 3, 'B': 1047552.0, 'pkdrop': 14}

Attack uses block-size 605 and 3072 dimensions, with 0 q-vectors
log2(epsilon) = -103.25, log2 nvector per run 125.55
shortest vector used has length l=24985036.35, q=8380417, `l<q'= 0
SIS & 3072 & 605 & 176 & 160 & 125
Primal attacks uses block-size 600 and 1195 samples
Primal & 1195 & 600 & 175 & 159 & 124
Dual attacks uses block-size 600 and 1150 samples
shortest vector used has length l=6191871.90, q=8380417, `l<q'= 1
log2(epsilon) = -62.18, log2 nvector per run 124.51
Dual & 1150 & 600 & 175 & 159 & 124
UNIFORM DILITHIUM TABLE
========================
\hline
BKZ block-size $b$ to break SIS &
235
 &
355
 &
475
 &
605
\\
Best Known Classical bit-cost &
68
 &
103
 &
138
 &
176
\\
Best Known Quantum bit-cost &
62
 &
94
 &
125
 &
160
\\
Best Plausible bit-cost &
48
 &
73
 &
98
 &
125
\\
\hline
BKZ block-size $b$ to break LWE &
200
 &
345
 &
485
 &
600
\\
Best Known Classical bit-cost &
58
 &
100
 &
141
 &
175
\\
Best Known Quantum bit-cost &
53
 &
91
 &
128
 &
159
\\
Best Plausible bit-cost &
41
 &
71
 &
100
 &
124
\\
\hline
========================
'''
