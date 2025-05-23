"""
Helper functions for the supersingular elliptic curve computations in FESTA
"""

# =========================================== #
# Compute points of order D and Torsion Bases #
# =========================================== #

"""
Most of this code has been taken from:
https://github.com/FESTA-PKE/FESTA-SageMath

Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope.

Functions with another Copyright mention are not from the above authors.
"""

# Sage Imports
from sage.all import ZZ, GF, Integers

# Local imports
from ..utilities.order import has_order_D
from ..utilities.discrete_log import weil_pairing_pari, ell_discrete_log_pari, discrete_log_pari
from ..utilities.fast_sqrt import sqrt_Fp2, sqrt_Fp




def random_point(E):
    """
    Returns a random point on the elliptic curve E
    assumed to be in Montgomery form with a base 
    field which characteristic p = 3 mod 4
    """
    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `generate_point` assumes the curve E is in the Montgomery model")

    # Try 10000 times then give up, just protection
    # for infinite loops
    F = E.base_ring()
    for _ in range(10000):
        x = F.random_element()
        y2 = x*(x**2 + A*x + 1)
        if y2.is_square():
            y = sqrt_Fp2(y2)
            return E(x, y)

    raise ValueError("Generated 10000 points, something is probably going wrong somewhere.")

def generate_point(E, x_start=0):
    """
    Generate points on a curve E with x-coordinate 
    i + x for x in Fp and i is the generator of Fp^2
    such that i^2 = -1.
    """
    F = E.base_ring()
    one = F.one()

    if x_start:
        x = x_start + one
    else:
        x = F.gen() + one

    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `generate_point` assumes the curve E is in the Montgomery model")

    # Try 10000 times then give up, just protection
    # for infinite loops
    for _ in range(10000):
        y2 = x*(x**2 + A*x + 1)
        if y2.is_square():
            y = sqrt_Fp2(y2)
            yield E(x, y)
        x += one

    raise ValueError("Generated 10000 points, something is probably going wrong somewhere.")

def generate_Fp_point(E):
    """
    Returns a random Fp-rational point on the 
    elliptic curve E assumed to be in Montgomery 
    form with a base field which characteristic 
    p = 3 mod 4 and coefficients defined over Fp.

    Copyright (c) Pierrick Dartois 2025.
    """

    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `generate_Fp_point` assumes the curve E is in the Montgomery model")
    if A[1] != 0:
        raise ValueError("The curve E should be Fp-rational")

    p=E.base_field().characteristic()
    Fp=GF(p,proof=False)
    one=Fp(1)
    x=one

    # Try 10000 times then give up, just protection
    # for infinite loops
    for _ in range(10000):
        y2 = x*(x**2 + A*x + 1)
        if y2.is_square():
            y = sqrt_Fp(y2)
            yield E(x, y)
        x += one

    raise ValueError("Generated 10000 points, something is probably going wrong somewhere.")

def generate_point_order_D(E, D, x_start=0):
    """
    Input:  An elliptic curve E / Fp2
            An integer D dividing (p +1)
    Output: A point P of order D.
    """
    p = E.base().characteristic()
    n = (p + 1) // D

    Ps = generate_point(E, x_start=x_start)
    for G in Ps:
        P = n * G

        # Case when we randomly picked
        # a point in the n-torsion
        if P.is_zero():
            continue

        # Check that P has order exactly D
        if has_order_D(P, D):
            P._order = ZZ(D)
            yield P

    raise ValueError(f"Never found a point P of order D.")

def generate_Fp_point_order_D(E, D):
    """
    Input:  An elliptic curve E / Fp2
            An integer D dividing (p +1)
    Output: A point P of order D.

    Copyright (c) Pierrick Dartois 2025.
    """
    p = E.base().characteristic()
    n = (p + 1) // D
    if D%2==0:
        n=n//2

    Ps = generate_Fp_point(E)
    for G in Ps:
        P = n * G

        # Case when we randomly picked
        # a point in the n-torsion
        if P.is_zero():
            continue

        # Check that P has order exactly D
        if has_order_D(P, D):
            P._order = ZZ(D)
            yield P

    raise ValueError(f"Never found a point P of order D.")

def compute_point_order_D(E, D, x_start=0):
    """
    Wrapper function around a generator which returns the first
    point of order D
    """
    return generate_point_order_D(E, D, x_start=x_start).__next__()

def compute_Fp_point_order_D(E, D):
    """
    Wrapper function around a generator which returns the first
    point of order D

    Copyright (c) Pierrick Dartois 2025.
    """
    return generate_Fp_point_order_D(E, D).__next__()

def compute_linearly_independent_point_with_pairing(E, P, D, x_start=0):
    """
    Input:  An elliptic curve E / Fp2
            A point P ∈ E[D]
            An integer D dividing (p +1)
    Output: A point Q such that E[D] = <P, Q>
            The Weil pairing e(P,Q)
    """
    Qs = generate_point_order_D(E, D, x_start=x_start)
    for Q in Qs:
        # Make sure the point is linearly independent
        pair = weil_pairing_pari(P, Q, D)
        if has_order_D(pair, D, multiplicative=True):
            Q._order = ZZ(D)
            return Q, pair
    raise ValueError("Never found a point Q linearly independent to P")

def compute_linearly_independent_point(E, P, D, x_start=0):
    """
    Wrapper function around `compute_linearly_independent_point_with_pairing`
    which only returns a linearly independent point
    """
    Q, _ = compute_linearly_independent_point_with_pairing(E, P, D, x_start=x_start)
    return Q

def torsion_basis_with_pairing(E, D):
    """
    Generate basis of E(Fp^2)[D] of supersingular curve

    While computing E[D] = <P, Q> we naturally compute the 
    Weil pairing e(P,Q), which we also return as in some cases
    the Weil pairing is then used when solving the BiDLP
    """
    p = E.base().characteristic()

    # Ensure D divides the curve's order
    if (p + 1) % D != 0:
        print(f"{ZZ(D).factor() = }")
        print(f"{ZZ(p+1).factor() = }")
        raise ValueError(f"D must divide the point's order")

    P = compute_point_order_D(E, D)
    Q, ePQ = compute_linearly_independent_point_with_pairing(E, P, D, x_start=P[0])

    return P, Q, ePQ

def torsion_basis(E, D):
    """
    Wrapper function around torsion_basis_with_pairing which only
    returns the torsion basis <P,Q> = E[D]
    """
    P, Q, _ = torsion_basis_with_pairing(E, D)
    return P, Q

def torsion_basis_2f_frobenius(E,f):
    """
    Returns a basis (P, Q) of E[2^f] with pi_p(P)=P and pi_p(Q)=-Q,
    pi_p being the p-th Frobenius endomorphism of E.

    Copyright (c) Pierrick Dartois 2025.
    """

    P = compute_Fp_point_order_D(E, 2**f)

    i = E.base_field().gen()
    p = E.base_field().characteristic() 

    R = compute_linearly_independent_point(E, P, 2**f)
    piR = E(R[0]**p,R[1]**p)
    piR_p_R = piR + R

    two_lamb = ell_discrete_log_pari(E,piR_p_R,P,2**f)
    lamb = two_lamb//2
    Q = R-lamb*P

    piQ = E(Q[0]**p,Q[1]**p)

    assert piQ == -Q

    return P, Q 

def torsion_basis_to_Fp_rational_point(E,P,Q,D):
    """
    Returns an Fp-rational point R=[lamb]P+[mu]Q
    of D-torsion given a D-torsion basis (P,Q) of E.

    Copyright (c) Pierrick Dartois 2025.
    """

    p=E.base_field().characteristic()
    piP = E(P[0]**p,P[1]**p)
    piQ = E(Q[0]**p,Q[1]**p)

    e_PQ = weil_pairing_pari(P,Q,D)
    e_PpiP = weil_pairing_pari(P, piP, D)
    e_QpiP = weil_pairing_pari(Q, piP, D)
    e_PpiQ = weil_pairing_pari(P, piQ, D)
    e_QpiQ = weil_pairing_pari(Q, piQ, D)

    # pi(P)=[a]P+[b]Q, pi(Q)=[c]P+[d]Q
    a = -discrete_log_pari(e_QpiP,e_PQ,D)
    b = discrete_log_pari(e_PpiP,e_PQ,D)
    c = -discrete_log_pari(e_QpiQ,e_PQ,D)
    d = discrete_log_pari(e_PpiQ,e_PQ,D)

    ZD = Integers(D)

    if a%D==1:
        mu = ZD(b/d)
        R=P+mu*Q
    elif c==0:
        # c=0 and a!=1 => d=1 so pi(Q)=Q
        R=Q
    else:
        # c!=0 and a!=1
        mu = ZD((1-a)/c)
        R=P+mu*Q

    piR = E(R[0]**p,R[1]**p)

    assert piR==R

    return R

# =========================================== #
#   Entangled torsion basis for fast E[2^k]   #
# =========================================== #

def precompute_elligator_tables(F):
    """
    Precomputes tables of quadratic non-residue or 
    quadratic residue in Fp2. Used to compute entangled
    torsion bases following https://ia.cr/2017/1143
    """
    u = 2*F.gen()

    T1 = dict()
    T2 = dict()
    # TODO: estimate how large r should be
    for r in range(1, 30):
        v = 1 / (1 + u*r**2)
        if v.is_square():
            T2[r] = v
        else:
            T1[r] = v
    return T1, T2

def entangled_torsion_basis(E, elligator_tables, cofactor):
    """
    Optimised algorithm following https://ia.cr/2017/1143
    which modifies the elligator method of hashing to points
    to find points P,Q of order k*2^b. Clearing the cofactor
    gives the torsion basis without checking the order or
    computing a Weil pairing. 

    To do this, we need tables TQNR, TQR of pairs of values
    (r, v) where r is an integer and v = 1/(1 + ur^2) where
    v is either a quadratic non-residue or quadratic residue 
    in Fp2 and u = 2i = 2*sqrt(-1).
    """
    F = E.base_ring()
    p = F.characteristic()
    p_sqrt = (p+1)//4

    i = F.gen()
    u =  2 * i
    u0 = 1 + i

    TQNR, TQR = elligator_tables
    
    # Pick the look up table depending on whether
    # A = a + ib is a QR or NQR
    A = E.a_invariants()[1]
    if (0,A,0,1,0) != E.a_invariants():
        raise ValueError("The elliptic curve E must be in Montgomery form")
    if A.is_square():
        T = TQNR
    else:
        T = TQR

    # Look through the table to find point with 
    # rational (x,y)
    y = None
    for r, v in T.items():
        x = -A * v

        t = x * (x**2 + A*x + 1)
        
        # Break when we find rational y: t = y^2 
        c, d = t.list()
        z = c**2 + d**2
        s = z**p_sqrt
        if s**2 == z:
            y = sqrt_Fp2(t)
            break
    
    if y is None:
        raise ValueError("Never found a y-coordinate, increase the lookup table size")

    z = (c + s) // 2
    alpha = z**p_sqrt
    beta  = d / (2*alpha)

    if alpha**2 == z:
        y =  F([alpha, beta])
    else:
        y = -F([beta, alpha])

    S1 = E([x, y])
    S2 = E([u*r**2*x, u0*r*y])

    return cofactor*S1, cofactor*S2

# =============================================== #
#  Ensure Basis <P,Q> of E[2^k] has (0,0) under Q #
# =============================================== #

def fix_torsion_basis_renes(P, Q, k):
    """
    Set the torsion basis P,Q such that
    2^(k-1)Q = (0,0) to ensure that (0,0)
    is never a kernel of a two isogeny
    """
    cofactor = 2**(k-1)

    R = cofactor*P
    if R[0] == 0:
        return Q, P
    R = cofactor*Q
    if R[0] == 0:
        return P, Q
    return P, P + Q



