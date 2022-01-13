from curve import BLS12
from proof import KZGProof
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element as vector

# for polynomials of degree 2^k
k = 2

# the number of accumulation
nb_accumulation = 2

# BLS12-381 curve
b = 4
alpha = -1
beta = [264,724]
pol = [593872, 0, 0, 0, 0, 0, 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559259, 0, 0, 0, 0, 0, 1]
z = -0xd201000000010000
Ebls = BLS12(z, pol, b, alpha, beta)         

####
# FFT setting

# ω_32 is a 2³²-th root of unity
ω_32 = Ebls.Fq(31368579106538355522267913380935626224263568375150176744936964454006709485244)
assert ω_32**(1<<32) == 1 and ω_32**(1<<31) != 1

# ω is a 2^n-th root of unity
n = 1<<k
ω = ω_32
for i in range(32-k):
    ω = ω**2
assert ω**n == 1 and ω**(1<<(k-1)) != 1

FX = Ebls.Fq['X']
X = FX.gen()
def lagrange_poly(i, size=n):
    pol = FX(1)
    for j in range(size):
        if j!= i:
            pol *= (X-ω**j) / (ω**i-ω**j)
    return pol

M =  Matrix(Ebls.Fq, n,n)
for i in range(n):
    Li = lagrange_poly(i).list()
    for j in range(n):
        M[j,i] = Li[j]
M_inv = M.inverse()

###

# s = Ebls.Fq.random_element()
s = Ebls.Fq(12)
trusted_setup = [int(s**i)*Ebls.G1 for i in range(1<<k)] + [Ebls.G2,int(s)*Ebls.G2]

P = KZGProof(trusted_setup, Ebls)
φ = Ebls.FqX([1,2,3,4])
# φ = Ebls.FqX([Ebls.Fq.random_element() for i in range(1<<k)])
P.commit(φ)

# a = Ebls.Fq.random_element()
a = Ebls.Fq(123)
y = φ(a)
P.open(φ, a)
assert P.verify()

#######################################

setup = [int(lagrange_poly(i, n)(s)) * Ebls.G1 for i in range(n)]

print('setup')
for point in setup:
    print(hex(point[0]))
    print(hex(point[1]))
    print()

φ_can = φ.list()
φ_lag = M_inv * vector(φ_can)
for coeff in φ_lag:
    print(hex(coeff))

P = KZGProof(setup+ [Ebls.G2,int(s)*Ebls.G2], Ebls)
P.commit(φ_lag)
print('commit')
print(hex(P.C[0]))
print(hex(P.C[1]))

y = sum([φ_lag[i] * lagrange_poly(i)(a) for i in range(n)])
P.open_lag(φ_lag, a, FX, M, M_inv)
print('proof')
print(hex(P.proof[0]))
print(hex(P.proof[1]))
assert P.verify()