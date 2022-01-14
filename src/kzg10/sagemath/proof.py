from copy import copy
from aux import rho
from sage.modules.free_module_element import free_module_element as vector

class KZGProof():
    def __init__(self, trusted_setup, curve):
        self.name='KZG'
        self.trusted_setup = trusted_setup
        self.C = 0
        self.proof = 0
        self.curve = curve
        self.a = 0
        self.y = 0

    def commit(self, φ, rand=0):
        # the commit is [φ(tau)]G, computed using the trusted setup.
        if isinstance(φ, list) :
            coeffs = φ
        else:
            coeffs = φ.list()
        assert len(coeffs) <= len(self.trusted_setup) - 2
        self.C = 0
        for i in range(len(coeffs)):
            self.C += int(coeffs[i]) * self.trusted_setup[i]
        return self.C

    def open(self, φ, a, rand=0):
        # computes a quotient polynomial and create the proof using
        # the trusted setup.
        y = φ(a)
        X = φ.parent().gen()
        q = (φ - y) // (X-a)
        coeffs = q.list()

        # assert len(coeffs) <= len(self.trusted_setup) - 2
        self.proof = 0
        for i in range(len(coeffs)):
            self.proof += int(coeffs[i]) * self.trusted_setup[i]
        self.a = a
        self.y = y
        return [self.y, self.proof]

    def open_lag(self, φ, a, ω):
        # computes a quotient polynomial and create the proof using
        # the trusted setup.

        # y = φ(a) is computed using the barycentric formula
        n = len(φ)
        res = 0
        pow_ω = 1
        for i in range(n):
            res += φ[i] * pow_ω / (pow_ω-a)
            pow_ω *= ω
        y = res * (1-a**n)/n
        coeffs = [(φ[i] - y)/(ω**i - a) for i in range(len(φ))]

        assert len(coeffs) <= len(self.trusted_setup) - 2
        self.proof = 0
        for i in range(len(coeffs)):
            self.proof += int(coeffs[i]) * self.trusted_setup[i]
        self.a = a
        self.y = y
        return [self.y, self.proof]
            
    def verify(self):
        C = self.curve.Ek(self.C)
        pi = self.curve.Ek(self.proof)
        G1 = self.curve.Ek(self.trusted_setup[-3])
        G2 = self.trusted_setup[-2]
        tauG2 = self.trusted_setup[-1]
        return (C - int(self.y)*G1).tate_pairing(G2, self.curve.q, 12)\
            == pi.tate_pairing(tauG2 - int(self.a)*G2, self.curve.q, 12)


