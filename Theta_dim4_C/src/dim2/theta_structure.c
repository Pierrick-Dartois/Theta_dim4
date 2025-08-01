#include "theta_structure.h"
#include <assert.h>

void
theta_precomputation(theta_structure_t *A)
{

    if (A->precomputation) {
        return;
    }

    theta_point_t A_dual;
    to_squared_theta(&A_dual, &A->null_point);

    fp2_t t1, t2;
    fp2_mul(&t1, &A_dual.x, &A_dual.y);
    fp2_mul(&t2, &A_dual.z, &A_dual.t);
    fp2_mul(&A->XYZ0, &t1, &A_dual.z);
    fp2_mul(&A->XYT0, &t1, &A_dual.t);
    fp2_mul(&A->YZT0, &t2, &A_dual.y);
    fp2_mul(&A->XZT0, &t2, &A_dual.x);

    fp2_mul(&t1, &A->null_point.x, &A->null_point.y);
    fp2_mul(&t2, &A->null_point.z, &A->null_point.t);
    fp2_mul(&A->xyz0, &t1, &A->null_point.z);
    fp2_mul(&A->xyt0, &t1, &A->null_point.t);
    fp2_mul(&A->yzt0, &t2, &A->null_point.y);
    fp2_mul(&A->xzt0, &t2, &A->null_point.x);

    A->precomputation = true;
}

static inline void
theta_semi_precomputation(theta_structure_t *A)
{

    if (A->semi_precomputation) {
        return;
    }

    theta_point_t A_dual;
    to_squared_theta(&A_dual, &A->null_point);

    fp2_t t1, t2;
    fp2_mul(&t1, &A_dual.x, &A_dual.y);
    fp2_mul(&t2, &A_dual.z, &A_dual.t);
    fp2_mul(&A->XYZ0, &t1, &A_dual.z);
    fp2_mul(&A->XYT0, &t1, &A_dual.t);
    fp2_mul(&A->YZT0, &t2, &A_dual.y);
    fp2_mul(&A->XZT0, &t2, &A_dual.x);

    A->semi_precomputation = true;
}

void
double_point(theta_point_t *out, theta_structure_t *A, const theta_point_t *in)
{
    to_squared_theta(out, in);
    fp2_sqr(&out->x, &out->x);
    fp2_sqr(&out->y, &out->y);
    fp2_sqr(&out->z, &out->z);
    fp2_sqr(&out->t, &out->t);

    if (!A->precomputation) {
        theta_precomputation(A);
    }
    fp2_mul(&out->x, &out->x, &A->YZT0);
    fp2_mul(&out->y, &out->y, &A->XZT0);
    fp2_mul(&out->z, &out->z, &A->XYT0);
    fp2_mul(&out->t, &out->t, &A->XYZ0);

    hadamard(out, out);

    fp2_mul(&out->x, &out->x, &A->yzt0);
    fp2_mul(&out->y, &out->y, &A->xzt0);
    fp2_mul(&out->z, &out->z, &A->xyt0);
    fp2_mul(&out->t, &out->t, &A->xyz0);
}

void
diff_add_point(theta_point_t *out, theta_structure_t *A, const theta_point_t *P, 
    const theta_point_t *Q, const theta_point_t *PQ)
{
    theta_point_t HSP, HSQ;

    to_squared_theta(&HSP, P);
    to_squared_theta(&HSQ, Q);

    field_mul(&out->x, &HSP.x, &HSQ.x);
    field_mul(&out->y, &HSP.y, &HSQ.y);
    field_mul(&out->z, &HSP.z, &HSQ.z);
    field_mul(&out->t, &HSP.t, &HSQ.t);

    if (!A->semi_precomputation) {// TODO: less precomp
        theta_semi_precomputation(A);
    }
    field_mul(&out->x, &out->x, &A->YZT0);
    field_mul(&out->y, &out->y, &A->XZT0);
    field_mul(&out->z, &out->z, &A->XYT0);
    field_mul(&out->t, &out->t, &A->XYZ0);

    hadamard(out, out);

    field_t xy, zt; 

    field_mul(&xy, &PQ->x, &PQ->y);
    field_mul(&zt, &PQ->z, &PQ->t);

    field_mul(&out->x, &out->x, &PQ->y);
    field_mul(&out->x, &out->x, &zt);

    field_mul(&out->y, &out->y, &PQ->x);
    field_mul(&out->y, &out->y, &zt);

    field_mul(&out->z, &out->z, &PQ->t);
    field_mul(&out->z, &out->z, &xy);

    field_mul(&out->t, &out->t, &PQ->z);
    field_mul(&out->t, &out->t, &xy);
}

void
double_iter(theta_point_t *out, theta_structure_t *A, const theta_point_t *in, int exp)
{
    if (exp == 0) {
        *out = *in;
    } else {
        double_point(out, A, in);
        for (int i = 1; i < exp; i++) {
            double_point(out, A, out);
        }
    }
}

uint32_t
is_product_theta_point(const theta_point_t *P)
{
    fp2_t t1, t2;
    fp2_mul(&t1, &P->x, &P->t);
    fp2_mul(&t2, &P->y, &P->z);
    return fp2_is_equal(&t1, &t2);
}
