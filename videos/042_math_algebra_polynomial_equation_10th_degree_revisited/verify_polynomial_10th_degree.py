"""Exact symbolic verification for the tenth-degree polynomial teaching note.

Run: python verify_polynomial_10th_degree.py

All arithmetic is exact (SymPy). Parameter sweeps test finite grids and do
not replace proofs of the general statements in the associated manuscript.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import gcd
import sympy as sp

x, t = sp.symbols("x t")
P = 1+x+x**2+x**4+x**5+x**6+x**8+x**9+x**10

@dataclass(frozen=True)
class CheckResult:
    name: str
    cases: int
    summary: str
    def __str__(self):
        return f"PASS  [{self.cases:4d} checks]  {self.name}: {self.summary}"

def _same(a, b):
    assert sp.cancel(a-b) == 0, f"Identity failed: {a} != {b}"

def _counts(f):
    entries = sp.Poly(f, x).sqf_list()[1]
    assert all(k in (1,2) for _,k in entries)
    simple = sum(g.degree() for g,k in entries if k==1)
    double = sum(g.degree() for g,k in entries if k==2)
    return simple+double, simple, double

def S(m, z=x):
    if m<1: raise ValueError("m must be positive")
    return sum((z**j for j in range(m)),sp.Integer(0))

def Q(m,n,b):
    return sp.expand(S(m)*S(n,x**b))

def Pmb(m,b):
    return Q(m,m,b)

def V(j,z=t):
    if j<0: raise ValueError("j must be nonnegative")
    a,b=sp.Integer(2),z
    if j==0: return a
    for _ in range(1,j): a,b=b,sp.expand(z*b-a)
    return b

def check_original_blocks():
    fact=S(3)*S(3,x**4)
    blocks=sum(x**(4*s)*S(3) for s in range(3))
    exps=[r+4*s for s in range(3) for r in range(3)]
    _same(P,fact); _same(P,blocks)
    assert len(exps)==len(set(exps))==9
    assert sorted(exps)==[0,1,2,4,5,6,8,9,10]
    assert sp.Poly(P,x).degree()==10
    return CheckResult("Exponent blocks",5,"support, degree and factorization")

def check_finite_geometric_sums():
    total=0
    for m in range(2,13):
        for b in range(1,9):
            z=x**b
            _same((z-1)*S(m,z),z**m-1)
            assert S(m,sp.Integer(1))==m
            total+=2
    return CheckResult("Finite geometric series",total,"telescoping and z=1")

def check_roots_of_unity():
    total=0
    for n in range(1,13):
        factors=sp.prod(sp.cyclotomic_poly(d,x) for d in sp.divisors(n))
        _same(x**n-1,factors)
        assert sum(sp.totient(d) for d in sp.divisors(n))==n
        for h in range(n):
            order=n//gcd(n,h)
            assert n%order==0
            assert (gcd(n,h)==1)==(order==n)
            total+=2
        total+=2
    return CheckResult("Roots of unity and orders",total,"order and divisor decomposition")

def check_cyclotomic_polynomials():
    expected={
        1:x-1,2:x+1,3:x*x+x+1,4:x*x+1,
        5:sum(x**j for j in range(5)),6:x*x-x+1,
        7:sum(x**j for j in range(7)),8:x**4+1,
        9:x**6+x**3+1,10:x**4-x**3+x*x-x+1,
        11:sum(x**j for j in range(11)),12:x**4-x*x+1,
    }
    for n,f in expected.items():
        phi=sp.Poly(sp.cyclotomic_poly(n,x),x)
        _same(phi.as_expr(),f)
        assert phi.degree()==sp.totient(n)
        assert phi.is_irreducible
    return CheckResult("Cyclotomic table Phi_1 through Phi_12",36,"forms, degrees, irreducibility")

def check_original_cyclotomic_factorization():
    p3,p6,p12=(sp.cyclotomic_poly(k,x) for k in (3,6,12))
    _same(P,p3**2*p6*p12)
    _same(1+x**4+x**8,p3*p6*p12)
    _same(1+x**4+x**8,(x*x+x+1)*(x*x-x+1)*(x**4-x*x+1))
    assert sp.Poly(sp.gcd(P,sp.diff(P,x)),x).as_expr()==p3
    return CheckResult("Original cyclotomic factorization",4,"Phi_3^2 Phi_6 Phi_12, gcd(P,P')")

def check_multiplicity_balance():
    _same((1-x)*(1-x**4)*P,(1-x**3)*(1-x**12))
    expected={h:(0 if h in (0,3,6,9) else 2 if h in (4,8) else 1) for h in range(12)}
    for h in range(12):
        balance=1+int(h%4==0)-int(h%3==0)-int(h==0)
        assert balance==expected[h]
    assert _counts(P)==(8,6,2)
    assert sum(expected.values())==10
    assert sp.Poly(sp.gcd(P,x**12-1),x).degree()==8
    return CheckResult("Twelve candidate roots and multiplicities",16,"four excluded, six simple, two double")

def check_fourth_power_substitution():
    u=sp.symbols("u")
    _same(1+x**4+x**8,(1+u+u**2).subs(u,x**4))
    plus=sorted((1+3*k)%12 for k in range(4))
    minus=sorted((-1+3*k)%12 for k in range(4))
    assert plus==[1,4,7,10]
    assert minus==[2,5,8,11]
    assert set(plus).isdisjoint(minus)
    assert {4,8}.issubset(set(plus+minus))
    return CheckResult("Substitution u=x^4",5,"eight fourth roots and two shared roots")

def check_real_factorization():
    a=x*x-sp.sqrt(3)*x+1
    b=x*x+sp.sqrt(3)*x+1
    _same(a*b,x**4-x*x+1)
    _same(P,(x*x+x+1)**2*(x*x-x+1)*a*b)
    for value in (-1,1,sp.sqrt(3),-sp.sqrt(3)):
        assert sp.simplify(value**2-4).is_negative
    return CheckResult("Real quadratic factorization",6,"all real quadratics have negative discriminant")

def check_viete():
    f=sp.Poly(P,x)
    assert [f.nth(j) for j in (10,9,8,0)]==[1,1,1,1]
    indices=[1,2,4,4,5,7,8,8,10,11]
    roots=[sp.exp(sp.I*sp.pi*sp.Rational(h,6)) for h in indices]
    sum1=sp.simplify(sp.expand_complex(sum(roots)))
    sum2=sp.simplify(sp.expand_complex(sum(r*r for r in roots)))
    product=sp.simplify(sp.expand_complex(sp.prod(roots)))
    assert sum1==-1 and sum2==-1 and product==1
    assert sp.simplify((sum1**2-sum2)/2)==1
    return CheckResult("Viete and power sums",5,"e1=-1, e2=1, e10=1")

def check_chebyshev_recurrences():
    y=sp.symbols("y")
    total=0
    for j in range(1,13):
        _same(sp.chebyshevt(j+1,y),2*y*sp.chebyshevt(j,y)-sp.chebyshevt(j-1,y))
        _same(sp.chebyshevu(j+1,y),2*y*sp.chebyshevu(j,y)-sp.chebyshevu(j-1,y))
        total+=2
    for j in range(11):
        _same(V(j),2*sp.chebyshevt(j,t/2))
        _same(x**j+x**(-j),V(j,x+1/x))
        total+=2
    return CheckResult("Chebyshev and symmetric powers",total,"T, U, V recurrences and Laurent identities")

def check_reciprocal_reduction():
    reduced=t**5+t**4-4*t**3-4*t*t+3*t+3
    expression=V(5)+V(4)+V(3)+V(1)+1
    factor=(t+1)**2*(t-1)*(t*t-3)
    _same(expression,reduced)
    _same(reduced,factor)
    _same(P,x**5*factor.subs(t,x+1/x))
    _same(P,x**5*expression.subs(t,x+1/x))
    coefficients=sp.Poly(P,x).all_coeffs()
    assert coefficients==list(reversed(coefficients))
    return CheckResult("Reciprocal quintic reduction",5,"fifth-degree reduction and its factorization")

def check_general_reciprocal_reduction():
    total=0
    for m in range(1,10):
        coefficients=[sp.Integer(j+1) for j in range(m)]
        f=x**m*(m+1)+sum(coefficients[j]*(x**j+x**(2*m-j)) for j in range(m))
        r=(m+1)+sum(coefficients[j]*V(m-j) for j in range(m))
        _same(f,x**m*r.subs(t,x+1/x))
        assert sp.Poly(r,t).degree()==m
        total+=2
    return CheckResult("General reciprocal reduction",total,"degree 2m to degree m")

def check_separated_exponent_blocks():
    total=0
    for m in range(2,8):
        for q in range(2,7):
            for b in range(m,m+5):
                exps=[r+b*s for s in range(q) for r in range(m)]
                assert len(set(exps))==m*q
                f=sum(x**e for e in exps)
                _same(f,S(m)*S(q,x**b))
                assert sp.Poly(f,x).degree()==m-1+b*(q-1)
                total+=3
    return CheckResult("Separated-block criterion",total,"support, factorization and degree")

def check_equal_length_counts():
    total=0
    for m in range(2,9):
        for b in range(1,9):
            f=Pmb(m,b)
            g=gcd(m,b)
            distinct,simple,double=_counts(f)
            assert sp.Poly(f,x).degree()==(m-1)*(b+1)
            assert double==m-g
            assert distinct==b*(m-1)+g-1
            assert simple+2*double==sp.Poly(f,x).degree()
            assert sp.Poly(sp.gcd(S(m),S(m,x**b)),x).degree()==m-g
            total+=5
    return CheckResult("Equal-length sums P_{m,b}",total,"gcd overlap and square-free root counts")

def check_unequal_length_counts():
    total=0
    for m in range(2,8):
        for n in range(2,8):
            for b in range(1,7):
                f=Q(m,n,b)
                delta=gcd(m,n*b)-gcd(m,b)
                distinct,simple,double=_counts(f)
                degree=(m-1)+b*(n-1)
                assert sp.Poly(f,x).degree()==degree
                assert double==delta
                assert distinct==degree-delta
                assert simple+2*double==degree
                assert sp.Poly(sp.gcd(S(m),S(n,x**b)),x).degree()==delta
                total+=5
    return CheckResult("Unequal-length sums Q_{m,n,b}",total,"gcd difference and root counts")

def check_real_root_parity():
    total=0
    for m in range(2,9):
        for n in range(2,9):
            for b in range(1,7):
                f=Q(m,n,b)
                multiplicity=int(m%2==0)+int(n%2==0 and b%2==1)
                assert (f.subs(x,-1)==0)==(multiplicity>0)
                if multiplicity:
                    assert sp.diff(f,x,multiplicity-1).subs(x,-1)==0
                    assert sp.diff(f,x,multiplicity).subs(x,-1)!=0
                else:
                    assert f.subs(x,0)==1 and f.subs(x,1)==m*n
                total+=3
    return CheckResult("Real-root parity criterion",total,"multiplicity of -1 is zero, one or two")

def check_worked_examples():
    examples=[
        ((3,3,2),(4,2,2,6)),((6,6,4),(21,17,4,25)),
        ((6,4,3),(11,8,3,14)),((4,4,2),(7,5,2,9)),
        ((3,3,3),(8,8,0,8)),((3,3,4),(8,6,2,10)),
        ((5,5,6),(24,20,4,28)),
    ]
    for (m,n,b),expected in examples:
        f=Q(m,n,b)
        assert (*_counts(f),sp.Poly(f,x).degree())==expected
    _same(Q(3,3,2),sp.cyclotomic_poly(3,x)**2*sp.cyclotomic_poly(6,x))
    _same(Q(4,4,2),(x+1)*(x*x+1)**2*(x**4+1))
    _same(Q(3,3,3),S(9))
    _same(sp.gcd(S(6),S(4,x**3)),(x+1)*(x*x-x+1))
    return CheckResult("Paper examples and four exercises",11,"seven examples and cyclotomic checks")

CHECKS=(
    check_original_blocks,check_finite_geometric_sums,check_roots_of_unity,
    check_cyclotomic_polynomials,check_original_cyclotomic_factorization,
    check_multiplicity_balance,check_fourth_power_substitution,
    check_real_factorization,check_viete,check_chebyshev_recurrences,
    check_reciprocal_reduction,check_general_reciprocal_reduction,
    check_separated_exponent_blocks,check_equal_length_counts,
    check_unequal_length_counts,check_real_root_parity,check_worked_examples,
)

def run_all(verbose=True):
    results=[check() for check in CHECKS]
    if verbose:
        for result in results: print(result)
        print(f"TOTAL: {sum(result.cases for result in results)} exact checks passed")
    return results

if __name__=="__main__":
    run_all()
