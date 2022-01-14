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
beta = [1,1]
pol = [2, 0,0,0,0,0, 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559785, 0,0,0,0,0,1]
z = -0xd201000000010000
BLS = BLS12(z, pol, b, alpha, beta)         

####
# FFT setting

# ω_32 is a 2³²-th root of unity
ω_32 = BLS.Fq(31368579106538355522267913380935626224263568375150176744936964454006709485244)
assert ω_32**(1<<32) == 1 and ω_32**(1<<31) != 1

# ω is a 2^n-th root of unity
n = 1<<k
ω = ω_32
for i in range(32-k):
    ω = ω**2
assert ω**n == 1 and ω**(1<<(k-1)) != 1

FX = BLS.Fq['X']
X = FX.gen()
def lagrange_poly(i, size=n):
    pol = FX(1)
    for j in range(size):
        if j!= i:
            pol *= (X-ω**j) / (ω**i-ω**j)
    return pol

M =  Matrix(BLS.Fq, n,n)
for i in range(n):
    Li = lagrange_poly(i).list()
    for j in range(n):
        M[j,i] = Li[j]
M_inv = M.inverse()

def evaluate_lag(φ_lag, a):
    # evaluate a polynomial φ given in Lagrange representation, at a using the barycentric formula
    n = len(φ_lag)
    res = 0
    pow_ω = BLS.Fq(1)
    for i in range(len(φ_lag)):
        res += φ_lag[i] * pow_ω / (pow_ω-a)
        pow_ω *= ω
    assert pow_ω == 1
    return res * (1-a**n)/n

###########################################
# these are the points obtained using the deterministic test_rng() in the rust code tests
# secret s
s = BLS.Fq(0x00E4C6AD736A4D69936C01D72DDF5848731A9C9249049CF6DD63FC02389BF946)
# G1
G1 = BLS.E(\
    0x07076CA73498ECAE77A45322438A3ED3633E63F8AE1A7F4FC22B9B830B16B2E9E6D2825382C4086465D869C2A1FD78F2,\
    0x0F0E2C6EF189E542C6DC0A906EF21445136EBC790E7C3A31BF49442E55D1FEBDD78112225E438EB54220DF8581588373\
    )
# G2
u = BLS.Fp2.gen()
G2 = BLS.E2_t(\
    0x194C1E61FB10CE68CE66C260C0E2A606B618784CF327E50FC30C3D5E98F5621F646A9077B8AD92DC59D4C007DF4B10C5 + \
        0x1244C9B4E93E42CB011A0C5D4EA369253A4608590F3B407006F452A7E4477A8CDE741B09980E06D96BE3540DDD95568C * u,
    0x0AF8AEC72D55D1DD3C0DDB6D14236DC65B706872AA42ADB330B38633F97E3829744F49C279188575004FA0605909B0CD + \
        0x12825EA037D0177B0B0ACD21F02CDB5AE433D1D719EC2234EDC0E4B964779C9FE27DB326B0BA5B137E86E259869E33F5 * u
    )
G2 = BLS.into_Ek(G2)
##########################################

# secret polynomial
φ = BLS.FqX([1,2,3,4])
# φ = BLS.FqX([BLS.Fq.random_element() for i in range(1<<k)])
φ_can = φ.list()
φ_lag = M_inv * vector(φ_can)

trusted_setup = [int(s**i)*G1 for i in range(1<<k)] + [G1, G2, int(s)*G2]
P1 = KZGProof(trusted_setup, BLS)
P1.commit(φ)
a1 = BLS.Fq(0x3F7EA3CC6271486D245360FD5B61B6E45D0E30CE3ABF2B77EBF2CB6D6FC8B94D) # as in test_rng() of arkworks
y1 = φ(a1)
P1.open(φ, a1)
assert P1.verify()

#######################################
trusted_setup_lagrange = [int(lagrange_poly(i, n)(s)) * G1 for i in range(n)] + [G1, G2,int(s)*G2]
P2 = KZGProof(trusted_setup_lagrange, BLS)
P2.commit(φ_lag)
a2 = a1
y2 = sum([φ_lag[i] * lagrange_poly(i)(a2) for i in range(n)])
P2.open_lag(φ_lag, a2, ω)
assert P2.verify()

assert P1.C == P2.C
assert P1.proof == P2.proof
