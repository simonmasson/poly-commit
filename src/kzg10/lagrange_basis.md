# Recall on Lagrange polynomials
Given our set {1, w, ..., w^{n-1}}, we define the j-th Lagrange polynomial: 
L_j(X) = \prod_{0<=m<n, m!=j} (X-w^m)/(w^j-w^m)

# Cost of the current trusted setup
Basically, the trusted setup is computed using `n` SM where each
element is `s` times the previous one

# Cost of the Lagrange basis trusted setup
Using the Lagrange polynomial basis, the trusted setup becomes:
{ L_0(s)G, ..., L_{n-1}(s)G }
In order to get this setup, we can:
1. Compute P = \prod_m (s-w^m)
2. For all j:
    * Divide P by (s-w^j) in order to get the right numerator,
    * Divide by prod_m (w^j-w^m)
	* Compute the scalar multiplication for G.
	
The overall complexity is roughly the same (if we look at scalar
multiplications).
