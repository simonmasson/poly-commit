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

    def commit(self, phi, rand=0):
        # the commit is [phi(tau)]G, computed using the trusted setup.
        if isinstance(phi, list) :
            coeffs = phi
        else:
            coeffs = phi.list()
        assert len(coeffs) <= len(self.trusted_setup) - 2
        self.C = 0
        for i in range(len(coeffs)):
            self.C += int(coeffs[i]) * self.trusted_setup[i]
        return self.C

    def open(self, phi, a, rand=0):
        # computes a quotient polynomial and create the proof using
        # the trusted setup.
        y = phi(a)
        X = phi.parent().gen()
        q = (phi - y) // (X-a)
        coeffs = q.list()   
        assert len(coeffs) <= len(self.trusted_setup) - 2
        self.proof = 0
        for i in range(len(coeffs)):
            self.proof += int(coeffs[i]) * self.trusted_setup[i]
        self.a = a
        self.y = y
        return [self.y, self.proof]

    def open_lag(self, phi, a, X, mat_inv):
        # computes a quotient polynomial and create the proof using
        # the trusted setup.

        # conversion in canonical basis
        phi_x = mat_inv.inverse() * vector(phi)
        # y = phi_x(a)
        y = sum([phi_x[i] * a**i for i in range(len(phi_x))])
        # q = (phi_x - y) / (X-a)
        q = (sum([phi_x[i] * X**i for i in range(len(phi))]) - y) // (X - a)
        vq = q.list()
        while len(vq)!= len(mat_inv.columns()):
            vq = vq + [0]
        coeffs = mat_inv * vector(vq)

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
        G1 = self.curve.Ek(self.curve.G1)
        G2 = self.trusted_setup[-2]
        tauG2 = self.trusted_setup[-1]
        return (C - int(self.y)*G1).tate_pairing(G2, self.curve.q, 12)\
            == pi.tate_pairing(tauG2 - int(self.a)*G2, self.curve.q, 12)
