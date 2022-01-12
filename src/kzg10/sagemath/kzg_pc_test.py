from curve import BLS12
from proof import KZGProof
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element as vector

# for polynomials of degree 2^k
k = 4

# the number of accumulation
nb_accumulation = 2

# BLS12-381 curve
b = 4
alpha = -1
beta = [264,724]
pol = [593872, 0, 0, 0, 0, 0, 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559259, 0, 0, 0, 0, 0, 1]
z = -0xd201000000010000
Ebls = BLS12(z, pol, b, alpha, beta)                                          

s = Ebls.Fq.random_element()
trusted_setup = [int(s**i)*Ebls.G1 for i in range(1<<k)] + [Ebls.G2,int(s)*Ebls.G2]

k = 4
P = KZGProof(trusted_setup, Ebls)
phi = Ebls.FqX([Ebls.Fq.random_element() for i in range(1<<k)])
P.commit(phi)
a = Ebls.Fq.random_element()
y = phi(a)
P.open(phi, a)
assert P.verify()

# ω_32 is a 2³²-th root of unity
ω_32 = Ebls.Fq(31368579106538355522267913380935626224263568375150176744936964454006709485244)
assert ω_32**(1<<32) == 1 and ω_32**(1<<31) != 1

# ω is a 2^n-th root of unity
n = 1<<k
ω = ω_32
for i in range(32-k):
    ω = ω**2
assert ω**n == 1 and ω**(1<<(k-1)) != 1

# TRUSTED SETUP IN LAGRANGE BASIS
num = Ebls.Fq(1)
den = Ebls.Fq(1)
pow_ω = 1
for i in range(n):
    num *= (s - pow_ω)
    if i>0:
	    den *= (1-pow_ω)
    pow_ω *= ω

setup = []
ω_n_minus_1 = ω**(n-1)
pow_ω = 1
for j in range(n):
    num_tmp = num / (s-pow_ω)
    setup.append(int(num_tmp/den) * Ebls.G1)
    pow_ω *= ω
    den *= ω_n_minus_1

# Check with a simple Lagrange polynomial construction
FX = Ebls.Fq['X']
X = FX.gen()
setup_check = []
def lagrange_poly(i):
    pol = FX(1)
    for j in range(n):
        if j!= i:
            pol *= (X-ω**j) / (ω**i-ω**j)
    return pol

M =  Matrix(Ebls.Fq, n,n)
for i in range(n):
    Li = lagrange_poly(i).list()
    for j in range(n):
        M[j,i] = Li[j]
M = M.inverse()


for i in range(n):
    poly = lagrange_poly(i)
    # check
    for j in range(n):
        if j!=i:
            assert poly(ω**j) == 0
        else :
            assert poly(ω**j) == 1

    setup_check.append(int(poly(s))*Ebls.G1)

assert setup == setup_check

phi_can = phi.list()
phi_lag = M * vector(phi_can)

k = 4
n = 1<<k
P = KZGProof(setup+ [Ebls.G2,int(s)*Ebls.G2], Ebls)
P.commit(phi_lag)
y = sum([phi_lag[i] * lagrange_poly(i)(a) for i in range(n)])
P.open_lag(phi_lag, a, X, M)

assert P.verify()
