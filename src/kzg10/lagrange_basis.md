# Recall on Lagrange polynomials

Given our set {1, ω, ..., ω^{n-1}}, we define the j-th Lagrange polynomial for the set
of n elements above as: 
L_{j, n}(X) = \prod_{0<=m<n, m!=j} (X-ω^m)/(ω^j-ω^m)

## Change of basis
The L_j form a basis of the polynomials of degree <= n and we can go from the
canonical representation to the Lagrangre representation using the matrix M
whose columns are the L_i in the canonical basis:
for a polynomial P = \sum_i a_i X^i, P_can = [a_i] and P_lag = [b_i] and
M * P_lag = P_can
M.inverse() * P_can = P_lag.
Note that P_lag can be obtained as [P(ω**i) for i in range(n)].

## Multiplication and division of polynomials
Given P(X) and Q(X) in Lagrange representation, say
P_lag = [P(ω**i) for i in range(n)]
Q_lag = [Q(ω**i) for i in range(n)]
the polynomial R(X) = P(X)*Q(X) has a Lagrange representation that
begins with
[P_lag[i] * Q_lag[i] for i in range(n)]
but its degree is greater than deg(P), deg(Q).

Given P(X) in Lagrange representation, and a degree 1 polynomial (X-a),
the representation of P(X)/(X-a) is
[P_lag[i] / (ω**i-a) for i in range(n-1)]
but because this is an interpolation with n-1 points, it corresponds to
the basis {1,ω,...,ω^{n-2}}.

## KZG in Lagrange basis
We modify the trusted setup so that we can apply the commitment scheme for 
polynomials given in Lagrange representation.
### Commitment
The commitment works exactly the same and there is no modification to be done.
### Opening proof
The proof opening requires computing the evaluation φ(a). Currently, the only way
to do that is to go back to the canonical representation of φ, in the base
1, X, ..., X^n.
The division polynomial could be done in the Lagrange representation. Though, in order
to "commit" in the last step of the opening proof, we  would need to have the 
{L_{i,n-1}(s) * G1, 0 <= i <= n-1} and it would increase a lot the size of the
trusted setup.  Moreover, we have already computed the canonical representation
for the evaluation φ(a), so we compute the division polynomial directly in the canonical
basis. In order to commit at the end of this step, we go back in Lagrange representation
using a fft.
In total, one ifft and one fft for the opening proof. It seems to be the same as if we do
it before the commitment scheme and work only with the canonical basis, but in the context
of PLONK, we have several polynomials and at the end only two opening proofs.
