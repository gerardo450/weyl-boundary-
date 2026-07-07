from __future__ import annotations
from dataclasses import dataclass
import sympy as sp
from math import factorial

@dataclass(frozen=True, order=True)
class Monomial:
    xdeg:int
    ddeg:int

class WeylOperator:
    def __init__(self,terms=None):
        self.terms={}
        if terms:
            for k,v in terms.items():
                if isinstance(k,tuple): k=Monomial(*k)
                v=sp.expand(v)
                if v!=0:self.terms[k]=v
    def copy(self): return WeylOperator(dict(self.terms))
    def simplify(self):
        self.terms={k:sp.expand(v) for k,v in self.terms.items() if sp.expand(v)!=0}; return self
    @staticmethod
    def _commute(m,n):
        return [(sp.binomial(m,k)*factorial(n)//factorial(n-k),n-k,m-k) for k in range(min(m,n)+1)]
    def __add__(self,o):
        if not isinstance(o,WeylOperator): o=o*I
        r=self.copy()
        for k,v in o.terms.items(): r.terms[k]=sp.expand(r.terms.get(k,0)+v)
        return r.simplify()
    __radd__=__add__
    def __neg__(self): return (-1)*self
    def __sub__(self,o): return self+(-o)
    def __rmul__(self,c): return WeylOperator({k:sp.expand(c*v) for k,v in self.terms.items()})
    def __mul__(self,o):
        if not isinstance(o,WeylOperator): return o*self
        out={}
        for m1,c1 in self.terms.items():
            for m2,c2 in o.terms.items():
                for cf,nx,nd in self._commute(m1.ddeg,m2.xdeg):
                    key=Monomial(m1.xdeg+nx,nd+m2.ddeg)
                    out[key]=sp.expand(out.get(key,0)+c1*c2*cf)
        return WeylOperator(out).simplify()
    def __pow__(self,n):
        r=I
        for _ in range(n): r=r*self
        return r
    def __eq__(self,o): return self.terms==o.terms
    def __repr__(self):
        items=sorted(self.terms.items(),key=lambda kv:(-(kv[0].xdeg+kv[0].ddeg),-kv[0].ddeg,-kv[0].xdeg))
        return " + ".join([str(c)+("" if m.xdeg==0 else ("*x" if m.xdeg==1 else f"*x^{m.xdeg}"))+("" if m.ddeg==0 else ("*D" if m.ddeg==1 else f"*D^{m.ddeg}")) for m,c in items]) or "0"
I=WeylOperator({Monomial(0,0):1})
X=WeylOperator({Monomial(1,0):1})
D=WeylOperator({Monomial(0,1):1})
def affine(a,b,c): return a*D+b*I+c*X

if __name__=="__main__":
    assert D*X==X*D+I
    assert D*(X**2)==(X**2)*D+2*X
    assert (D**2)*X==X*(D**2)+2*D
    print("tests passed")
