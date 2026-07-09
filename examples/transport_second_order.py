import sympy as sp

import weyl.transport as wt
from weyl.algebra import *

# ----------------------------------------------------
# Symbols
# ----------------------------------------------------

t = sp.Symbol("t")

a0, a1, a2 = sp.symbols("a0 a1 a2")
b0, b1     = sp.symbols("b0 b1")
c0          = sp.symbols("c0")

# ----------------------------------------------------
# Build the operator
# ----------------------------------------------------

L = (
    c0*(D**2)
    +
    (b0*I + b1*X)*D
    +
    a0*I
    +
    a1*X
    +
    a2*(X**2)
)

print("Original operator")
print(L)
print()

Lt = wt.heat_transport(L, t)

print("Transported operator")
for m, c in sorted(Lt.terms.items()):
    print(m, ":", c)