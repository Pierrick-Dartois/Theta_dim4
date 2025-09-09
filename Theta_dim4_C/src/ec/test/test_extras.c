#include "test_extras.h"
#include "rng.h"

// Given an x-coordinate, determines if this is a valid
// point on the curve. Assumes C=1.
static uint32_t
projective_is_on_curve(const ec_point_t *P, const ec_curve_t *curve)
{

    field_t t0, t1, t2;

    // Check if xz*(C^2x^2+zACx+z^2C^2) is a square
    field_mul(&t0, &curve->C, &P->x);
    field_mul(&t1, &t0, &P->z);
    field_mul(&t1, &t1, &curve->A);
    field_mul(&t2, &curve->C, &P->z);
    field_sqr(&t0, &t0);
    field_sqr(&t2, &t2);
    field_add(&t0, &t0, &t1);
    field_add(&t0, &t0, &t2);
    field_mul(&t0, &t0, &P->x);
    field_mul(&t0, &t0, &P->z);
    return field_is_square(&t0) || field_is_zero(&t0);
}

void
ec_random_normalized_test(ec_point_t *P, const ec_curve_t *curve)
{
    field_set_one(&P->z);
    while (1) {
        field_random_test(&P->x);
        if (projective_is_on_curve(P, curve)) {
            break;
        }
    }
}

void
ec_random_test(ec_point_t *P, const ec_curve_t *curve)
{
    ec_random_normalized_test(P, curve);
    field_random_test(&P->z);
    field_mul(&P->x, &P->x, &P->z);
}

