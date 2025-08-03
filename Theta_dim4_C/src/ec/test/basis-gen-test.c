#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <ec.h>
#include "test_extras.h"

/******************************
Test functions
******************************/

int
inner_test_generated_basis(ec_basis_t *basis, ec_curve_t *curve, unsigned int n)
{
    unsigned int i;
    int passed = 1;

    ec_point_t P, Q;
    copy_point(&P, &basis->P);
    copy_point(&Q, &basis->Q);

    // Double points to get point of order 2
    for (i = 0; i < n - 1; i++) {
        xDBL_A24(&P, &P, &curve->A24, curve->is_A24_computed_and_normalized);
        xDBL_A24(&Q, &Q, &curve->A24, curve->is_A24_computed_and_normalized);
    }
    if (ec_is_zero(&P)) {
        printf("Point P generated does not have full order\n");
        passed = 0;
    }
    if (ec_is_zero(&Q)) {
        printf("Point Q generated does not have full order\n");
        passed = 0;
    }
    if (ec_is_equal(&P, &Q)) {
        printf("Points P, Q are linearly dependent\n");
        passed = 0;
    }

    if (!fp2_is_zero(&Q.x)) {
        printf("Points Q is not above the Montgomery point\n");
        passed = 0;
    }

    // This should give the identity
    xDBL_A24(&P, &P, &curve->A24, curve->is_A24_computed_and_normalized);
    xDBL_A24(&Q, &Q, &curve->A24, curve->is_A24_computed_and_normalized);
    if (!ec_is_zero(&P)) {
        printf("Point P generated does not have order exactly 2^n\n");
        passed = 0;
    }
    if (!ec_is_zero(&Q)) {
        printf("Point Q generated does not have order exactly 2^n\n");
        passed = 0;
    }

    if (passed == 0) {
        printf("Test failed with n = %u\n", n);
    }

    return passed;
}

int
inner_test_hint_basis(ec_basis_t *basis, ec_basis_t *basis_hint)
{
    int passed = 1;

    if (!ec_is_equal(&basis->P, &basis_hint->P)) {
        printf("The points P do not match using the hint\n");
        passed = 0;
    }

    if (!ec_is_equal(&basis->Q, &basis_hint->Q)) {
        printf("The points Q do not match using the hint\n");
        passed = 0;
    }

    if (!ec_is_equal(&basis->PmQ, &basis_hint->PmQ)) {
        printf("The points PmQ do not match using the hint\n");
        passed = 0;
    }

    if (passed == 0) {
        printf("Test failed\n");
    }

    return passed;
}

/******************************
Test wrapper functions
******************************/

int
test_basis_generation_E0(unsigned int n)
{
    ec_basis_t basis;
    ec_curve_t curve;

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 0);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    (void)ec_curve_to_basis_2f_to_hint(&basis, &curve, n);

    // Test result
    return inner_test_generated_basis(&basis, &curve, n);
}

int
test_basis_generation(unsigned int n)
{
    ec_basis_t basis;
    ec_curve_t curve;

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    (void)ec_curve_to_basis_2f_to_hint(&basis, &curve, n);

    // Test result
    return inner_test_generated_basis(&basis, &curve, n);
}

int
test_basis_generation_with_hints(unsigned int n)
{
    int check_1, check_2;
    ec_basis_t basis, basis_hint;
    ec_curve_t curve;
    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis with hints
    uint8_t hint = ec_curve_to_basis_2f_to_hint(&basis, &curve, n);

    // Ensure the basis from the hint is good
    check_1 = inner_test_generated_basis(&basis, &curve, n);

    // Generate a basis using hints
    ec_curve_to_basis_2f_from_hint(&basis_hint, &curve, n, hint);

    // These two bases should be the same
    check_2 = inner_test_hint_basis(&basis, &basis_hint);

    return check_1 && check_2;
}

int
test_point_diff_invariance(int reps)
{
    ec_point_t P, Q, PQ, PQ2;
    ec_curve_t curve;

    for(int i = 0; i < reps; i++) {
        fp2_random_test(&curve.A);
        fp2_random_test(&curve.C);
        ec_random_test(&P, &curve);
        ec_random_test(&Q, &curve);

        projective_difference_point(&PQ, &P, &Q, &curve);

        ec_normalize_curve(&curve);
        ec_normalize_point(&P);
        ec_normalize_point(&Q);

        projective_difference_point(&PQ2, &P, &Q, &curve);

        fp2_mul(&PQ.x, &PQ.x, &PQ2.z);
        fp2_mul(&PQ2.x, &PQ2.x, &PQ.z);

        if(!fp2_is_equal(&PQ.x, &PQ2.x)) {
            return 0;
        }
    }
    return 1;
}

int
test_basis(void)
{
    int passed;

    // Test full order
    passed = test_basis_generation(TORSION_EVEN_POWER);
    passed &= test_basis_generation_with_hints(TORSION_EVEN_POWER);

    // Test partial order
    passed &= test_basis_generation(128);
    passed &= test_basis_generation_with_hints(128);

    // Special case when we have A = 0
    passed &= test_basis_generation_E0(TORSION_EVEN_POWER);
    passed &= test_basis_generation_E0(128);

    // Test that point difference does not depend on projective representation
    passed &= test_point_diff_invariance(1024);

    return passed;
}

int
main(void)
{
    bool ok;
    ok = test_basis();
    if (!ok) {
        printf("Tests failed!\n");
    } else {
        printf("All basis generation tests passed.\n");
    }
    return !ok;
}
