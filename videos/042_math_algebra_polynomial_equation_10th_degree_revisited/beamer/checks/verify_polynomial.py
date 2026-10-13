"""Exact symbolic checks for the 15-slide Beamer presentation.

Run: python checks/verify_polynomial.py
"""
import sympy as sp

x = sp.symbols('x')
P = 1 + x + x**2 + x**4 + x**5 + x**6 + x**8 + x**9 + x**10
Q = (1 + x + x**2)*(1 + x**4 + x**8)

checks = {
    'elementary factorization': sp.expand(P - Q) == 0,
    'geometric-series identity': sp.expand((1-x)*(1-x**4)*P - (1-x**3)*(1-x**12)) == 0,
    'cyclotomic factorization': sp.expand(P - sp.cyclotomic_poly(3,x)**2*sp.cyclotomic_poly(6,x)*sp.cyclotomic_poly(12,x)) == 0,
    'factor degrees and multiplicities': sp.Poly(P).degree() == 10 and sorted((sp.degree(f), power) for f,power in sp.factor_list(P)[1]) == [(2,1),(2,2),(4,1)],
    'rational factorization': sp.expand(P - (x**2+x+1)**2*(x**2-x+1)*(x**4-x**2+1)) == 0,
    'conjugate quadratic merging': sp.expand((x**2-sp.sqrt(3)*x+1)*(x**2+sp.sqrt(3)*x+1)-(x**4-x**2+1)) == 0,
    'sum of multiplicities': sum(sp.degree(f)*power for f,power in sp.factor_list(P)[1]) == 10,
    'no real roots': sp.Poly(P,x).count_roots(-sp.oo, sp.oo) == 0,
}
for name, ok in checks.items():
    print(f'{name}: {"PASS" if ok else "FAIL"}')
if not all(checks.values()):
    raise SystemExit('Some checks failed')
