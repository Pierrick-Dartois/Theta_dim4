
# ============================================ #
#     Fast square root and quadratic roots     #
# ============================================ #

"""
Most of this code has been taken from:
https://github.com/FESTA-PKE/FESTA-SageMath

Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope.

Functions with another Copyright mention are not from the above authors.
"""

def sqrt_Fp2(a):
    """
    Efficiently computes the sqrt
    of an element in Fp2 using that
    we always have a prime p such that
    p ≡ 3 mod 4.
    """
    Fp2 = a.parent()
    p = Fp2.characteristic()
    i = Fp2.gen() # i = √-1

    a1 = a ** ((p - 3) // 4)
    x0 = a1 * a
    alpha = a1 * x0

    if alpha == -1:
        x = i * x0
    else:
        b = (1 + alpha) ** ((p - 1) // 2)
        x = b * x0

    return x

def n_sqrt(a, n):
    for _ in range(n):
        a = sqrt_Fp2(a)
    return a

def sqrt_Fp(a):
    """
    Efficiently computes the sqrt
    of an element in Fp using that
    we always have a prime p such that
    p ≡ 3 mod 4.

    Copyright (c) Pierrick Dartois 2025.
    """
    Fp = a.parent()
    p = Fp.characteristic()

    return a**((p+1)//4)

# Deterministic square roots
def sqrt_Fp_det(a):
    """
    Efficiently computes the sqrt
    of an element in Fp using that
    we always have a prime p such that
    p ≡ 3 mod 4 and does it in a 
    deterministic way.

    Copyright (c) Pierrick Dartois 2025.
    """
    Fp = a.parent()
    p = Fp.characteristic()
    t = a**((p+1)//4)

    # Negate if necessary to make sure 
    # the sqrt is always even.
    if int(t)%2:
        t=-t

    return t

def sqrt_Fp2_det(a):
    """
    Efficiently computes the sqrt
    of an element in Fp2 using that
    we always have a prime p such that
    p ≡ 3 mod 4 and does it in a 
    deterministic way.

    Translation of sqrt_Fp2 function from 
    SQIsign2D-West repository.

    Copyright (c) Pierrick Dartois 2025.
    """
    Fp2 = a.parent()
    p = Fp2.characteristic()
    i = Fp2.gen() # i = √-1


    sqrt_delta = sqrt_Fp_det(a[0]**2 + a[1]**2)
    if a[1]==0:
        y02 = a[0]
    else:
        y02 = (a[0]+sqrt_delta)/2

    nqr = (y02**((p-1)//2)==-1)# 1 iff y02 is not a square in Fp

    if nqr:
        if a[1]==0:
            y02 = -y02
        else:
            # Take the other root (a[0]-sqrt_delta)/2
            y02 = y02-sqrt_delta

    y0 = sqrt_Fp_det(y02)
    y1 = a[1]/(2*y0)

    # If x1 = 0 then the sqrt worked and y1 = 0, but if
    # y02 was not a square, we must swap y0 and y1
    if a[0]==0 and nqr:
        y0, y1 = y1, y0

    # To ensure deterministic square-root we conditionally negate the output
    # We negate the output if:
    # y0 is odd, or
    # y0 is zero and y1 is odd
    if int(y0)%2 or (y0==0 and int(y1)%2):
        y = -y0-i*y1
    else:
        y = y0+i*y1

    return y