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

# s = Ebls.Fq.random_element()
s = Ebls.Fq(12)
trusted_setup = [int(s**i)*Ebls.G1 for i in range(1<<k)] + [Ebls.G2,int(s)*Ebls.G2]

P = KZGProof(trusted_setup, Ebls)
phi = Ebls.FqX([1,2,3,4])
# phi = Ebls.FqX([Ebls.Fq.random_element() for i in range(1<<k)])
P.commit(phi)
a = Ebls.Fq.random_element()
y = phi(a)
P.open(phi, a)

assert P.verify()

#######################################

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
# lagrange polynomial evaluations are num/den_i/(s-ω**i) where:
#  * num = prod_i (s-ω**i)
#  * den_i = prod_i (1-ω**i) * (ω**(n-1))
# Hence, we compute num, ω_n_minus_1, and den = prod_i(1-ω**i) first.
num = Ebls.Fq(1)
den = Ebls.Fq(1)
pow_ω = 1
for i in range(n):
    num *= (s - pow_ω)
    if i>0:
	    den *= (1-pow_ω)
    if i == n-1:
        ω_n_minus_1 = pow_ω
    pow_ω *= ω

setup = []
pow_ω = 1
for j in range(n):
    num_tmp = num / (s-pow_ω)
    setup.append(int(num_tmp/den) * Ebls.G1)
    pow_ω *= ω
    den *= ω_n_minus_1

# Check with a simple Lagrange polynomial construction
FX = Ebls.Fq['X']
X = FX.gen()
def lagrange_poly(i, size=n):
    pol = FX(1)
    for j in range(size):
        if j!= i:
            pol *= (X-ω**j) / (ω**i-ω**j)
    return pol

# could be more efficient as we did for `setup`
setup_small = [int(lagrange_poly(i,n-1)(s)) * Ebls.G1 for i in range(n-1)]

M =  Matrix(Ebls.Fq, n,n)
for i in range(n):
    Li = lagrange_poly(i).list()
    for j in range(n):
        M[j,i] = Li[j]
M = M.inverse()

for i in range(n):
    assert setup[i] == int(lagrange_poly(i)(s))*Ebls.G1

phi_can = phi.list()
phi_lag = M * vector(phi_can)

P = KZGProof(setup+setup_small+ [Ebls.G2,int(s)*Ebls.G2], Ebls)
P.commit(phi_lag)
y = sum([phi_lag[i] * lagrange_poly(i)(a) for i in range(n)])
P.open_lag(phi_lag, a, \
    [lagrange_poly(i) for i in range(n)],\
        [lagrange_poly(i,size=3) for i in range(3)],\
            X, ω, M)

assert P.verify()

x1 = 0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB
y1 = 0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1
x2 = 0x0345DD80FFEF0EAEC8920E39EBB7F5E9AE9C1D6179E9129B705923DF7830C67F3690CBC48649D4079EADF5397339580C
y2 = 0x083D3BAF25E42F2845D8FA594DDA2E0F40A4D670DDA40F30DA0AFF0D81C87AC3D687FE84ECA72F34C7C755A045668CF1
x3 = 0x07DC2DA68D1641FFE8E6CA1B675767DC3303995C5E9E31564905C196E3109F11345B8877D28D116E8AE110E6A6A7C7A4
y3 = 0x00C0BE655B2AA556FEAE51F9F43668B9FE66821A1F6A18E4C3E0C0D27FF4E6B29A0A5A7CFBB6D101BDBA2EF24070B495

x11 = 0x16A2E13765C1771BE1D19F7C41A38A82877163E31BB6BE20A341091081E16682C7874442950068EB7767A4E445289CAC
y11 = 0x07C55ADAFD29FFFDCD2327D35178659D51D28C01EC893602F823845C955F160CCEE1ADAC9934B476DDCE3E3CE5ABF9AA
x21 = 0x118CC9B07A9B42847FE5D12E185E72DC6B6AB6A4A1F893FE88EBAEB24A9104D34EC67AC93E96C7E9FDA8FD1C8C11CC06
y21 = 0x100CC6E2CEDB935517972B797BF4D72B765687624B3FE2543AFCC913AFA45B65C78A0E31BBD0B1AE277C0C26AB96C637
x31 = 0x07941E7B359221F84729D78A3F5F9F1CDA1894F6BBFF74683E945256A630B80F7D27D51E061A76149352BBBABFA6D66A
y31 = 0x02F32295344B637302C3597DD64A3B5C59F30F4DD4DB8D74F59A0491A4F2525134374CA7AA9A89EF3F71F696F0A7896C

assert trusted_setup[:3] == [Ebls.E(x1,y1), Ebls.E(x2,y2), Ebls.E(x3,y3)]
assert setup[:3] == [Ebls.E(x11,y11), Ebls.E(x21,y21), Ebls.E(x31,y31)]