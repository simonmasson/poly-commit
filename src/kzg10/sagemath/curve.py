# -*- coding: utf-8 -*-
from sage.all import ZZ
from sage.rings.rational_field import Q as QQ
# from sage.rings.integer_ring import Z as ZZ
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.functions.other import sqrt
from sage.rings.number_field.number_field import NumberField
from sage.structure.proof.all import arithmetic as ProofArithmetic

ProofArithmetic(False)

class Degree12PairingCurve():
    def __init__(self, p, q, k, PolFpk, alpha=False, beta=False, b=False):
        # alpha is a non-square of Fp defining Fp2. PolFpk defines Fp12/Fp
        self.p = p
        self.q = q
        self.Fp = GF(p)
        self.Fq = GF(q)
        self.FqX = self.Fq['X']
        self.k = k

        if b:
            self.E = EllipticCurve([self.Fp(0), self.Fp(b)])
        if not(b):
            b = self.Fp(1)
            self.E = EllipticCurve([self.Fp(0), b])
            while self.E.order() % self.q != 0 :
                b += 1
                self.E = EllipticCurve([self.Fp(0), b])
                self.b = b

        FpT = self.Fp['T']
        T = FpT.gen()
        if not(alpha):
            alpha = self.Fp.random_element()
            polFp2 = T**2 - alpha
            while not(polFp2.is_irreducible()):
                alpha = self.Fp.random_element()
                polFp2 = T**2 - alpha
        alpha = self.Fp(alpha)
        self.polFp2 = T**2 - alpha

        self.Fp2 = GF(self.p**2, 'u', modulus = self.polFp2)
        self.u = self.Fp2.gen()
        
        if not(beta):
            Fp2Z = self.Fp2['Z']
            Z = Fp2Z.gen()
            beta = self.Fp2.random_element()
            polFq6 = Z**6 - beta
            while not(polFq6.is_irreducible()):
                beta = self.Fp2.random_element()
                polFq6 = Z**6 - beta
        self.beta = self.Fp2(beta)

        self.G1 = self.E(3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507, 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569)
        # orderE = self.E.order()
        # cof = self.E.order()//self.q
        # self.G1 = cof*self.E.random_point()
        # while self.G1 == 0:
        #     self.G1 = cof*self.E.random_point()
        # assert self.G1!= 0 and self.q*self.G1==0

        self.E2 = self.E.change_ring(self.Fp2)

        self.E2_t = EllipticCurve([self.Fp2(0), self.beta * self.b])
        self.twist_elt = self.beta
        if self.E2_t.order() % self.q != 0:
            self.twist_elt = self.twist_elt**5
            self.E2_t = EllipticCurve([self.Fp2(0), self.twist_elt * self.b])
        assert self.E2_t.order() % self.q == 0

        cof2 = self.E2_t.order() // self.q
        self.G2 = cof2 * self.E2_t.random_point()
        while self.G2 == 0:
            self.G2 = cof2 * self.E2_t.random_point()

        self.Fpk = GF(self.p**12, 'v', modulus=PolFpk)
        self.v = self.Fpk.gen()
        assert self.beta.minpoly() == (self.v**6).minpoly()

        # The sextic twist is defined using beta:
        # E_{0,b} --> E_{0, twist_elt * b}
        # (x,y)  |--> (twist_elt.sqrt3()*x, twist_elt.sqrt()*y)
        Elt = self.into_Fpk(self.twist_elt)
        self.sqrt_elt = Elt.sqrt()
        Fpkz = self.Fpk['zz']
        zz = Fpkz.gen()
        self.sqrt3_elt = (zz**3 - Elt).roots()[0][0]
        self.Ek = self.E.change_ring(self.Fpk)

        self.G2 = self.into_Ek(self.G2)
        
    def into_Fpk(self, elt):
        # return an Fp2 element in the Fp12 field (which is a
        # Fp-extension)
        # Fp -- x²+1 -- Fp2 -- x⁶-beta -- Fp12
        a = self.beta.polynomial().list()
        if len(a) == 1 :
            a = a + [0]
        e = elt.polynomial().list()
        if len(e) == 1:
            e = e + [0]
        return self.Fp(e[0]) + self.Fp(e[1]) * (self.v**6-self.Fp(a[0])) / self.Fp(a[1])

    def into_Ek(self, pt):
            # return the twisted point to a point of E2_t (i.e. in Ek)
            x = self.into_Fpk(pt[0])
            y = self.into_Fpk(pt[1])
            return self.Ek([x / self.sqrt3_elt, y / self.sqrt_elt])



class BLS12(Degree12PairingCurve):
    def __init__(self, z, polFp12, b=False, alpha=False, beta=False):
        '''
        alpha is a non-square of Fp, defining Fp2
        beta is a non 6-th root of Fp2, defining Fp12 over Fp2
        polFp12 defines Fp12 over Fp, with a well defined injection
        Fp2 --> Fp12.
        '''
        C = ZZ['x']
        x = C.gen()
        p_x = (x-1)**2 * (x**4-x**2+1) + 3*x # needs a division by 3!
        q_x = (x**4-x**2+1)

        super().__init__(p_x(z)//3, q_x(z), 12, polFp12, alpha, beta) 
