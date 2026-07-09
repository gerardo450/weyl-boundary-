"""
transport.py

Heat transport in the Weyl algebra.

Implements

    X -> X + t D
    D -> D
"""

import sympy as sp

from .algebra import WeylOperator, X, D, I


def substitute(operator, Ximage):
    """
    Replace

        X -> Ximage
        D -> D

    in a Weyl operator.

    Parameters
    ----------
    operator : WeylOperator

    Ximage : WeylOperator
        Usually X + t*D.

    Returns
    -------
    WeylOperator
    """

    result = 0*I

    for monomial, coeff in operator.terms.items():

        piece = coeff * (Ximage**monomial.xdeg) * (D**monomial.ddeg)

        result = result + piece

    return result


def heat_transport(operator, t=None):
    """
    Compute

        T_t L T_t^{-1}

    where

        X -> X+tD.
    """

    if t is None:
        t = sp.Symbol("t")

    Xt = X + t*D

    return substitute(operator, Xt)


if __name__ == "__main__":

    t = sp.Symbol("t")

    from .algebra import affine

    L = affine(2,3,4)

    print("Original")
    print(L)

    print()

    print("Transported")
    print(heat_transport(L,t))