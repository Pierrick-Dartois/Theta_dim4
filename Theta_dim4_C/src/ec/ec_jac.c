#include <assert.h>
#include <ec.h>

void
jac_init(jac_point_t *P)
{ // Initialize Montgomery in Jacobian coordinates as identity element (0:1:0)
    field_set_zero(&P->x);
    field_set_one(&P->y);
    field_set_zero(&P->z);
}

uint32_t
jac_is_equal(const jac_point_t *P, const jac_point_t *Q)
{ // Evaluate if two points in Jacobian coordinates (X:Y:Z) are equal
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    field_t t0, t1, t2, t3;

    field_sqr(&t0, &Q->z);
    field_mul(&t2, &P->x, &t0); // x1*z2^2
    field_sqr(&t1, &P->z);
    field_mul(&t3, &Q->x, &t1); // x2*z1^2
    field_sub(&t2, &t2, &t3);

    field_mul(&t0, &t0, &Q->z);
    field_mul(&t0, &P->y, &t0); // y1*z2^3
    field_mul(&t1, &t1, &P->z);
    field_mul(&t1, &Q->y, &t1); // y2*z1^3
    field_sub(&t0, &t0, &t1);

    return field_is_zero(&t0) & field_is_zero(&t2);
}

void
jac_to_xz(ec_point_t *P, const jac_point_t *xyP)
{
    field_copy(&P->x, &xyP->x);
    field_copy(&P->z, &xyP->z);
    field_sqr(&P->z, &P->z);

    // If xyP = (0:1:0), we currently have P=(0 : 0) but we want to set P=(1:0)
    uint32_t c1, c2;
    field_t one;
    field_set_one(&one);

    c1 = field_is_zero(&P->x);
    c2 = field_is_zero(&P->z);
    field_select(&P->x, &P->x, &one, c1 & c2);
}

void
jac_to_ws(jac_point_t *Q, field_t *t, field_t *ao3, const jac_point_t *P, const ec_curve_t *curve)
{
    // Cost of 3M + 2S when A != 0.
    field_t one;
    field_t a;
    /* a = 1 - A^2/3, U = X + (A*Z^2)/3, V = Y, W = Z, T = a*Z^4*/
    field_set_one(&one);
    if (!field_is_zero(&(curve->A))) {
        field_div3(&ao3, &curve->A);
        field_sqr(t, &P->z);
        field_mul(&Q->x, ao3, t);
        field_add(&Q->x, &Q->x, &P->x);
        field_sqr(t, t);
        field_mul(&a, ao3, &(curve->A));
        field_sub(&a, &one, &a);
        field_mul(t, t, &a);
    } else {
        field_copy(&Q->x, &P->x);
        field_sqr(t, &P->z);
        field_sqr(t, t);
    }
    field_copy(&Q->y, &P->y);
    field_copy(&Q->z, &P->z);
}

void
jac_from_ws(jac_point_t *Q, const jac_point_t *P, const field_t *ao3, const ec_curve_t *curve)
{
    // Cost of 1M + 1S when A != 0.
    field_t t;
    /* X = U - (A*W^2)/3, Y = V, Z = W. */
    if (!field_is_zero(&(curve->A))) {
        field_sqr(&t, &P->z);
        field_mul(&t, &t, ao3);
        field_sub(&Q->x, &P->x, &t);
    }
    field_copy(&Q->y, &P->y);
    field_copy(&Q->z, &P->z);
}

void
copy_jac_point(jac_point_t *P, const jac_point_t *Q)
{
    field_copy(&(P->x), &(Q->x));
    field_copy(&(P->y), &(Q->y));
    field_copy(&(P->z), &(Q->z));
}

void
jac_neg(jac_point_t *Q, const jac_point_t *P)
{
    field_copy(&Q->x, &P->x);
    field_neg(&Q->y, &P->y);
    field_copy(&Q->z, &P->z);
}

void
DBL(jac_point_t *Q, const jac_point_t *P, const ec_curve_t *AC)
{ // Cost of 6M + 6S.
  // Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
  // (X/Z^2,Y/Z^3) This version receives the coefficient value A
    field_t t0, t1, t2, t3;

    uint32_t flag = field_is_zero(&P->x) & field_is_zero(&P->z);

    field_sqr(&t0, &P->x); // t0 = x1^2
    field_add(&t1, &t0, &t0);
    field_add(&t0, &t0, &t1); // t0 = 3x1^2
    field_sqr(&t1, &P->z);    // t1 = z1^2
    field_mul(&t2, &P->x, &AC->A);
    field_add(&t2, &t2, &t2); // t2 = 2Ax1
    field_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    field_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    field_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)
    field_mul(&Q->z, &P->y, &P->z);
    field_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1
    field_sqr(&t0, &Q->z);
    field_mul(&t0, &t0, &AC->A); // t0 = 4Ay1^2z1^2
    field_sqr(&t1, &P->y);
    field_add(&t1, &t1, &t1);     // t1 = 2y1^2
    field_add(&t3, &P->x, &P->x); // t3 = 2x1
    field_mul(&t3, &t1, &t3);     // t3 = 4x1y1^2
    field_sqr(&Q->x, &t2);        // x2 = alpha^2
    field_sub(&Q->x, &Q->x, &t0); // x2 = alpha^2 - 4Ay1^2z1^2
    field_sub(&Q->x, &Q->x, &t3);
    field_sub(&Q->x, &Q->x, &t3); // x2 = alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    field_sub(&Q->y, &t3, &Q->x); // y2 = 4x1y1^2 - x2
    field_mul(&Q->y, &Q->y, &t2); // y2 = alpha(4x1y1^2 - x2)
    field_sqr(&t1, &t1);          // t1 = 4y1^4
    field_sub(&Q->y, &Q->y, &t1);
    field_sub(&Q->y, &Q->y, &t1); // y2 = alpha(4x1y1^2 - x2) - 8y1^4

    field_select(&Q->x, &Q->x, &P->x, -flag);
    field_select(&Q->z, &Q->z, &P->z, -flag);
}

void
DBLW(jac_point_t *Q, field_t *u, const jac_point_t *P, const field_t *t)
{ // Cost of 3M + 5S.
  // Doubling on a Weierstrass curve, representation in modified Jacobian coordinates
  // (X:Y:Z:T=a*Z^4) corresponding to (X/Z^2,Y/Z^3), where a is the curve coefficient.
  // Formula from https://hyperelliptic.org/EFD/g1p/auto-shortw-modified.html

    uint32_t flag = field_is_zero(&P->x) & field_is_zero(&P->z);

    field_t xx, c, cc, r, s, m;
    // XX = X^2
    field_sqr(&xx, &P->x);
    // A = 2*Y^2
    field_sqr(&c, &P->y);
    field_add(&c, &c, &c);
    // AA = A^2
    field_sqr(&cc, &c);
    // R = 2*AA
    field_add(&r, &cc, &cc);
    // S = (X+A)^2-XX-AA
    field_add(&s, &P->x, &c);
    field_sqr(&s, &s);
    field_sub(&s, &s, &xx);
    field_sub(&s, &s, &cc);
    // M = 3*XX+T1
    field_add(&m, &xx, &xx);
    field_add(&m, &m, &xx);
    field_add(&m, &m, t);
    // X3 = M^2-2*S
    field_sqr(&Q->x, &m);
    field_sub(&Q->x, &Q->x, &s);
    field_sub(&Q->x, &Q->x, &s);
    // Z3 = 2*Y*Z
    field_mul(&Q->z, &P->y, &P->z);
    field_add(&Q->z, &Q->z, &Q->z);
    // Y3 = M*(S-X3)-R
    field_sub(&Q->y, &s, &Q->x);
    field_mul(&Q->y, &Q->y, &m);
    field_sub(&Q->y, &Q->y, &r);
    // T3 = 2*R*T1
    field_mul(u, t, &r);
    field_add(u, u, u);

    field_select(&Q->x, &Q->x, &P->x, -flag);
    field_select(&Q->z, &Q->z, &P->z, -flag);
}

void
select_jac_point(jac_point_t *Q, const jac_point_t *P1, const jac_point_t *P2, const digit_t option)
{ // Select points
  // If option = 0 then Q <- P1, else if option = 0xFF...FF then Q <- P2
    field_select(&(Q->x), &(P1->x), &(P2->x), option);
    field_select(&(Q->y), &(P1->y), &(P2->y), option);
    field_select(&(Q->z), &(P1->z), &(P2->z), option);
}

void
ADD(jac_point_t *R, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *AC)
{
    // Addition on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding
    // to (x,y) = (X/Z^2,Y/Z^3) This version receives the coefficient value A
    //
    // Complete routine, to handle all edge cases:
    //   if ZP == 0:            # P == inf
    //       return Q
    //   if ZQ == 0:            # Q == inf
    //       return P
    //   dy <- YQ*ZP**3 - YP*ZQ**3
    //   dx <- XQ*ZP**2 - XP*ZQ**2
    //   if dx == 0:             # x1 == x2
    //       if dy == 0:         # ... and y1 == y2: doubling case
    //           dy <- ZP*ZQ * (3*XP^2 + ZP^2 * (2*A*XP + ZP^2))
    //           dx <- 2*YP*ZP
    //       else:              # ... but y1 != y2, thus P = -Q
    //           return inf
    //   XR <- dy**2 - dx**2 * (A*ZP^2*ZQ^2 + XP*ZQ^2 + XQ*ZP^2)
    //   YR <- dy * (XP*ZQ^2 * dx^2 - XR) - YP*ZQ^3 * dx^3
    //   ZR <- dx * ZP * ZQ

    // Constant time processing:
    // - The case for P == 0 or Q == 0 is handled at the end with conditional select
    // - dy and dx are computed for both the normal and doubling cases, we switch when
    //   dx == dy == 0 for the normal case.
    // - If we have that P = -Q then dx = 0 and so ZR will be zero, giving us the point
    //   at infinity for "free".
    //
    // These current formula are expensive and I'm probably missing some tricks...
    // Thought I'd get the ball rolling.
    // Cost 17M + 6S + 13a
    field_t t0, t1, t2, t3, u1, u2, v1, dx, dy;

    /* If P is zero or Q is zero we will conditionally swap before returning. */
    uint32_t ctl1 = field_is_zero(&P->z);
    uint32_t ctl2 = field_is_zero(&Q->z);

    /* Precompute some values */
    field_sqr(&t0, &P->z); // t0 = z1^2
    field_sqr(&t1, &Q->z); // t1 = z2^2

    /* Compute dy and dx for ordinary case */
    field_mul(&v1, &t1, &Q->z); // v1 = z2^3
    field_mul(&t2, &t0, &P->z); // t2 = z1^3
    field_mul(&v1, &v1, &P->y); // v1 = y1z2^3
    field_mul(&t2, &t2, &Q->y); // t2 = y2z1^3
    field_sub(&dy, &t2, &v1);   // dy = y2z1^3 - y1z2^3
    field_mul(&u2, &t0, &Q->x); // u2 = x2z1^2
    field_mul(&u1, &t1, &P->x); // u1 = x1z2^2
    field_sub(&dx, &u2, &u1);   // dx = x2z1^2 - x1z2^2

    /* Compute dy and dx for doubling case */
    field_add(&t1, &P->y, &P->y);   // dx_dbl = t1 = 2y1
    field_add(&t2, &AC->A, &AC->A); // t2 = 2A
    field_mul(&t2, &t2, &P->x);     // t2 = 2Ax1
    field_add(&t2, &t2, &t0);       // t2 = 2Ax1 + z1^2
    field_mul(&t2, &t2, &t0);       // t2 = z1^2 * (2Ax1 + z1^2)
    field_sqr(&t0, &P->x);          // t0 = x1^2
    field_add(&t2, &t2, &t0);       // t2 = x1^2 + z1^2 * (2Ax1 + z1^2)
    field_add(&t2, &t2, &t0);       // t2 = 2*x1^2 + z1^2 * (2Ax1 + z1^2)
    field_add(&t2, &t2, &t0);       // t2 = 3*x1^2 + z1^2 * (2Ax1 + z1^2)
    field_mul(&t2, &t2, &Q->z);     // dy_dbl = t2 = z2 * (3*x1^2 + z1^2 * (2Ax1 + z1^2))

    /* If dx is zero and dy is zero swap with double variables */
    uint32_t ctl = field_is_zero(&dx) & field_is_zero(&dy);
    field_select(&dx, &dx, &t1, ctl);
    field_select(&dy, &dy, &t2, ctl);

    /* Some more precomputations */
    field_mul(&t0, &P->z, &Q->z); // t0 = z1z2
    field_sqr(&t1, &t0);          // t1 = z1z2^2
    field_sqr(&t2, &dx);          // t2 = dx^2
    field_sqr(&t3, &dy);          // t3 = dy^2

    /* Compute x3 = dy**2 - dx**2 * (A*ZP^2*ZQ^2 + XP*ZQ^2 + XQ*ZP^2) */
    field_mul(&R->x, &AC->A, &t1); // x3 = A*(z1z2)^2
    field_add(&R->x, &R->x, &u1);  // x3 = A*(z1z2)^2 + u1
    field_add(&R->x, &R->x, &u2);  // x3 = A*(z1z2)^2 + u1 + u2
    field_mul(&R->x, &R->x, &t2);  // x3 = dx^2 * (A*(z1z2)^2 + u1 + u2)
    field_sub(&R->x, &t3, &R->x);  // x3 = dy^2 - dx^2 * (A*(z1z2)^2 + u1 + u2)

    /* Compute y3 = dy * (XP*ZQ^2 * dx^2 - XR) - YP*ZQ^3 * dx^3*/
    field_mul(&R->y, &u1, &t2);     // y3 = u1 * dx^2
    field_sub(&R->y, &R->y, &R->x); // y3 = u1 * dx^2 - x3
    field_mul(&R->y, &R->y, &dy);   // y3 = dy * (u1 * dx^2 - x3)
    field_mul(&t3, &t2, &dx);       // t3 = dx^3
    field_mul(&t3, &t3, &v1);       // t3 = v1 * dx^3
    field_sub(&R->y, &R->y, &t3);   // y3 = dy * (u1 * dx^2 - x3) - v1 * dx^3

    /* Compute z3 = dx * z1 * z2 */
    field_mul(&R->z, &dx, &t0);

    /* Finally, we need to set R = P is Q.Z = 0 and R = Q if P.Z = 0 */
    select_jac_point(R, R, Q, ctl1);
    select_jac_point(R, R, P, ctl2);
}

void
jac_to_xz_add_components(add_components_t *add_comp, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *AC)
{
    // Take P and Q in E distinct, two jac_point_t, return three components u,v and w in the base field such
    // that the xz coordinates of P+Q are (u-v:w) and of P-Q are (u+v:w)

    field_t t0, t1, t2, t3, t4, t5, t6;

    field_sqr(&t0, &P->z);             // t0 = z1^2
    field_sqr(&t1, &Q->z);             // t1 = z2^2
    field_mul(&t2, &P->x, &t1);        // t2 = x1z2^2
    field_mul(&t3, &t0, &Q->x);        // t3 = z1^2x2
    field_mul(&t4, &P->y, &Q->z);      // t4 = y1z2
    field_mul(&t4, &t4, &t1);          // t4 = y1z2^3
    field_mul(&t5, &P->z, &Q->y);      // t5 = z1y2
    field_mul(&t5, &t5, &t0);          // t5 = z1^3y2
    field_mul(&t0, &t0, &t1);          // t0 = (z1z2)^2
    field_mul(&t6, &t4, &t5);          // t6 = (z1z_2)^3y1y2
    field_add(&add_comp->v, &t6, &t6); // v  = 2(z1z_2)^3y1y2
    field_sqr(&t4, &t4);               // t4 = y1^2z2^6
    field_sqr(&t5, &t5);               // t5 = z1^6y_2^2
    field_add(&t4, &t4, &t5);          // t4 = z1^6y_2^2 + y1^2z2^6
    field_add(&t5, &t2, &t3);          // t5 = x1z2^2 +z_1^2x2
    field_add(&t6, &t3, &t3);          // t6 = 2z_1^2x2
    field_sub(&t6, &t5, &t6);          // t6 = lambda = x1z2^2 - z_1^2x2
    field_sqr(&t6, &t6);               // t6 = lambda^2 = (x1z2^2 - z_1^2x2)^2
    field_mul(&t1, &AC->A, &t0);       // t1 = A*(z1z2)^2
    field_add(&t1, &t5, &t1);          // t1 = gamma =A*(z1z2)^2 + x1z2^2 +z_1^2x2
    field_mul(&t1, &t1, &t6);          // t1 = gamma*lambda^2
    field_sub(&add_comp->u, &t4, &t1); // u  = z1^6y_2^2 + y1^2z2^6 - gamma*lambda^2
    field_mul(&add_comp->w, &t6, &t0); // w  = (z1z2)^2(lambda)^2
}
