from sage.all import *
from ..montgomery_isogenies.kummer_line import KummerLine
from ..utilities.fast_sqrt import sqrt_Fp2_det


def difference_point(KP,KQ,E):
    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `difference_point` assumes the curve E is in the Montgomery model")

    xP = KP.x()
    xQ = KQ.x()

    zPQ = xP-xQ # P-Q
    t2 = xP*xQ # P*Q
    t3 = t2-1 # P*Q-1
    t0 = zPQ*t3 # (P-Q)*(P*Q-1)
    zPQ = zPQ*zPQ # (P-Q)^2
    t0 = t0*t0 # (P-Q)^2*(P*Q-1)^2
    t1 = t2+1 # P*Q+1
    t3 = xP+xQ # P+Q
    t1 = t1*t3 # (P+Q)*(P*Q+1)
    t2 = t2*A # A*P*Q
    t2 = t2 + t2 # 2*A*P*Q
    t1 = t1 + t2 # (P+Q)*(P*Q+1) + 2*A*P*Q
    t2 = t1*t1 # ((P+Q)*(P*Q+1) + 2*A*P*Q)^2
    t0 = t2-t0 # ((P+Q)*(P*Q+1) + 2*A*P*Q)^2 - (P-Q)^2*(P*Q-1)^2
    t0 = sqrt_Fp2_det(t0)
    xPQ = t0 + t1

    K = KP.parent()
    KPQ = K((xPQ,zPQ))

    return KPQ

def recover_y(x,E):
    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `recover_y` assumes the curve E is in the Montgomery model")

    y2 = x*(x*x+A*x+1)
    return sqrt_Fp2_det(y2)

def lift_basis(KP,KQ,KPQ,E):
    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `lift_basis` assumes the curve E is in the Montgomery model")

    xP = KP.x()
    xQ = KQ._X
    zQ = KQ._Z
    xPQ = KPQ._X
    zPQ = KPQ._Z 

    yP = recover_y(xP,E)

    # Algorithm of Okeya-Sakurai to recover y.Q in the montgomery model
    v1 = xP*zQ # fp2_mul(&v1, &P->x, &Q->z);
    v2 = xQ+v1 # fp2_add(&v2, &Q->x, &v1);
    v3 = xQ-v1 # fp2_sub(&v3, &Q->x, &v1);
    v3 = v3*v3 # fp2_sqr(&v3, &v3);
    v3 = v3*xPQ # fp2_mul(&v3, &v3, &B->PmQ.x);
    v1 = A+A # fp2_add(&v1, &E->A, &E->A);
    v1 = v1*zQ # fp2_mul(&v1, &v1, &Q->z);
    v2 = v2+v1 # fp2_add(&v2, &v2, &v1);
    v4 = xP*xQ # fp2_mul(&v4, &P->x, &Q->x);
    v4 = v4+zQ # fp2_add(&v4, &v4, &Q->z);
    v2 = v2*v4 # fp2_mul(&v2, &v2, &v4);
    v1 = v1*zQ # fp2_mul(&v1, &v1, &Q->z);
    v2 = v2-v1 # fp2_sub(&v2, &v2, &v1);
    v2 = v2*zPQ # fp2_mul(&v2, &v2, &B->PmQ.z);
    yQ = v3-v2 # fp2_sub(&Q->y, &v3, &v2);
    v1 = yP+yP # fp2_add(&v1, &P->y, &P->y);
    v1 = v1*zQ # fp2_mul(&v1, &v1, &Q->z);
    v1 = v1*zPQ # fp2_mul(&v1, &v1, &B->PmQ.z);
    xQ = xQ*v1 # fp2_mul(&Q->x, &Q->x, &v1);
    zQ = zQ*v1 # fp2_mul(&Q->z, &Q->z, &v1);

    return (E([xP,yP]), E([xQ,yQ,zQ]))

# ======================================== #
#  Generate a 2^f-torsion basis from hint  #
# ======================================== #

def torsion_basis_2f_from_hint(E,hP,hQ,NQR_TABLE,Z_NQR_TABLE):
    A = E.a_invariants()[1]
    if E.a_invariants() != (0,A,0,1,0):
        raise ValueError("Function `torsion_basis_2f_from_hint` assumes the curve E is in the Montgomery model")

    Fp2 = E.base_field()
    p = Fp2.characteristic()
    i = Fp2.gen()
    f = valuation(p+1,2)
    c = (p+1)//2**f
    K = KummerLine(E)

    if hP<20:
        xP = NQR_TABLE[hP]
    else:
        xP = hP + i

    alpha = (-A+sqrt_Fp2_det(A**2-4))/2
    if hQ<20:
        xQ = Z_NQR_TABLE[hQ]*alpha
    else:
        xQ = (hQ + i)*alpha
    KP = K(xP)
    KQ = K(xQ)

    KP = c*KP
    KQ = c*KQ

    KPQ = difference_point(KP,KQ,E)

    xP = KP.x()
    xQ = KP.x()
    K0 = K.zero()

    A = E.a_invariants()[1]

    return lift_basis(KP,KQ,KPQ,E)




    


