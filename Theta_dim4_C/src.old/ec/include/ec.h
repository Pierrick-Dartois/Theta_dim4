﻿/** @file
 *
 * @authors Luca De Feo, Francisco RH
 *
 * @brief Elliptic curve stuff
*/

#ifndef EC_H
#define EC_H

#include "../../gf/include/fp2.h"
//#include <ec_params.h>


/** @defgroup ec Elliptic curves
 * @{
*/

/** @defgroup ec_t Data structures
 * @{
*/

/** @brief Projective point
 *
 * @typedef ec_point_t
 *
 * @struct ec_point_t
 *
 * A projective point in (X:Z) or (X:Y:Z) coordinates (tbd).
*/
typedef struct ec_point_t {
    fp2_t x;
    fp2_t z;
} ec_point_t;

/** @brief A basis of a torsion subgroup
 *
 * @typedef ec_basis_t
 *
 * @struct ec_basis_t
 *
 * A pair of points (or a triplet, tbd) forming a basis of a torsion subgroup.
*/
typedef struct ec_basis_t {
    ec_point_t P;
    ec_point_t Q;
    ec_point_t PmQ;
} ec_basis_t;

/** @brief An elliptic curve
 *
 * @typedef ec_curve_t
 *
 * @struct ec_curve_t
 *
 * An elliptic curve in projective Montgomery form
*/
typedef struct ec_curve_t {
    fp2_t A;
    fp2_t C; ///< cannot be 0
} ec_curve_t;

/** @brief Projective Jacobian point
 *
 * @typedef jac_point_t
 *
 * @struct jac_point_t
 *
 * A projective point (X:Y:Z) in Jaconian form.
*/
typedef struct jac_point_t {
    fp2_t x;
    fp2_t y;
    fp2_t z;
} jac_point_t;

/** @brief An isogeny of degree a power of 2
 *
 * @typedef ec_isog_even_t
 *
 * @struct ec_isog_even_t
*/
typedef struct ec_isog_even_t {
    ec_curve_t curve;      ///< The domain curve
    ec_point_t kernel;     ///< A kernel generator
    unsigned short length; ///< The length as a 2-isogeny walk
} ec_isog_even_t;


/** @brief An odd divisor of p² - 1
 *
 * @typedef ec_isog_odd_t
 *
 * Given that the list of divisors of p² - 1 is known, this is
 * represented as a fixed-length vector of integer exponents.
*/

//typedef uint8_t ec_degree_odd_t[P_LEN + M_LEN];

/** @brief An isogeny of odd degree dividing p² - 1
 *
 * @typedef ec_isog_odd_t
 *
 * @struct ec_isog_odd_t
*/
//typedef struct ec_isog_odd_t {
    //ec_curve_t curve;
    //ec_point_t ker_plus;    ///< A generator of E[p+1] ∩ ker(φ)
    //ec_point_t ker_minus;   ///< A generator of E[p-1] ∩ ker(φ)
    //ec_degree_odd_t degree; ///< The degree of the isogeny
//} ec_isog_odd_t;

/** @brief Isomorphism of Montgomery curves
 *
 * @typedef ec_isom_t
 *
 * @struct ec_isom_t
 *
 * The isomorphism is given by the map maps (X:Z) ↦ ( (Nx X - Nz Z) : (D Z) )
*/
typedef struct ec_isom_t {
    fp2_t Nx;
    fp2_t Nz;
    fp2_t D;
} ec_isom_t;

// end ec_t
/** @}
*/


/** @defgroup ec_curve_t Curves and isomorphisms
 * @{
*/

/**
 * @brief j-invariant.
 *
 * @param j_inv computed j_invariant
 * @param curve input curve
 */
void ec_j_inv(fp2_t* j_inv, const ec_curve_t* curve);

/**
 * @brief Isomorphism of elliptic curve
 *
 * @param isom computed isomorphism
 * @param from domain curve
 * @param to image curve
 */
void ec_isomorphism(ec_isom_t* isom, const ec_curve_t* from, const ec_curve_t* to);

/**
 * @brief In-place inversion of an isomorphism
 *
 * @param isom an isomorphism
 */
void ec_iso_inv(ec_isom_t* isom);

/**
 * @brief In-place evaluation of an isomorphism
 *
 * @param P a point
 * @param isom an isomorphism
 */
void ec_iso_eval(ec_point_t* P, ec_isom_t* isom);

/**
 * @brief Given a Montgomery curve, computes a standard model for it and the isomorphism to it.
 *
 * @param new computed new curve
 * @param isom computed isomorphism from `old` to `new`
 * @param old A Montgomery curve
 */
void ec_curve_normalize(ec_curve_t *new, ec_isom_t *isom, const ec_curve_t *old);

/** @}
*/
/** @defgroup ec_point_t Point operations
 * @{
*/

/**
 * @brief Point equality
 *
 * @param P a point
 * @param Q a point
 * @return 1 if equal
 */
bool ec_is_equal(const ec_point_t* P, const ec_point_t* Q);

/**
 * @brief Reduce Z-coordinate of point in place
 *
 * @param P a point
 */
void ec_normalize(ec_point_t* P);

/**
 * @brief Test whether a point is on a curve
 *
 * @param curve a curve
 * @param P a point
 * @return 1 if P is on the curve
 */
int ec_is_on_curve(const ec_curve_t* curve, const ec_point_t* P);

/**
 * @brief Point negation
 *
 * @param res computed opposite of P
 * @param P a point
 */
void ec_neg(ec_point_t* res, const ec_point_t* P);

/**
 * @brief Point addition
 *
 * @param res computed sum of P and Q
 * @param P a point
 * @param Q a point
 * @param PQ the difference P-Q
 */
void ec_add(ec_point_t* res, const ec_point_t* P, const ec_point_t* Q, const ec_point_t* PQ);

/**
 * @brief Point doubling
 *
 * @param res computed double of P
 * @param P a point
 */
void ec_dbl(ec_point_t* res, const ec_curve_t* curve, const ec_point_t* P);

/**
 * @brief Point multiplication
 *
 * @param res computed scalar * P
 * @param curve the curve
 * @param scalar an unsigned multi-precision integer
 * @param P a point
 */
void ec_mul(ec_point_t* res, const ec_curve_t* curve, const digit_t* scalar, const ec_point_t* P);

/**
 * @brief Point multiplication by a scalar of limited length
 *
 * @param res computed scalar * P
 * @param curve the curve
 * @param scalar an unsigned multi-precision integer
 * @param kbits the bit size of scalar
 * @param P a point
 */
void xMULv2(ec_point_t* Q, ec_point_t const* P, digit_t const* k, const int kbits, ec_point_t const* A24);

/**
 * @brief Combination P+m*Q
 *
 * @param R computed P + m * Q
 * @param curve the curve
 * @param m an unsigned multi-precision integer
 * @param P a point
 * @param Q a point
 * @param PQ the difference P-Q
 */
void ec_ladder3pt(ec_point_t *R, digit_t const* m, ec_point_t const *P, ec_point_t const *Q, ec_point_t const *PQ, ec_curve_t const *A);

/**
 * @brief Linear combination of points of a basis
 *
 * @param res computed scalarP * P + scalarQ * Q
 * @param curve the curve
 * @param scalarP an unsigned multi-precision integer
 * @param scalarQ an unsigned multi-precision integer
 * @param PQ a torsion basis consisting of points P and Q
 */
void ec_biscalar_mul(ec_point_t* res, const ec_curve_t* curve,
    const digit_t* scalarP, const digit_t* scalarQ,
    const ec_basis_t* PQ);

/** @}
*/

/** @defgroup ec_dlog_t Discrete logs and bases
 * @{
*/

/**
 * @brief Generate a Montgomery curve and a 2^f-torsion basis
 *
 * The algorithm is deterministc
 *
 * @param PQ2 computed basis of the 2^f-torsion
 * @param curve the computed curve
 */
void ec_curve_to_basis_2(ec_basis_t *PQ2, const ec_curve_t *curve);

/**
 * @brief Complete a basis of the 2^f-torsion
 *
 * The algorithm is deterministic
 *
 * @param PQ2 a basis of the 2^f-torsion containing P as first generator
 * @param curve the curve
 * @param P a point of order 2^f
 */
void ec_complete_basis_2(ec_basis_t* PQ2, const ec_curve_t* curve, const ec_point_t* P);

/**
 * @brief Generate a 3^e-torsion basis
 *
 * The algorithm is deterministic
 *
 * @param PQ3 the computed 3^e-torsion basis
 * @param curve a curve
 */
void ec_curve_to_basis_3(ec_basis_t* PQ3, const ec_curve_t* curve);

/**
 * @brief Generate a 6^e-torsion basis
 *
 * The algorithm is deterministic
 *
 * @param PQ6 the computed 2^f*3^g-torsion basis
 * @param curve a curve
 */
void ec_curve_to_basis_6(ec_basis_t* PQ6, const ec_curve_t* curve);


/**
 * @brief Compute the generalized dlog of R wrt the 2^f-basis PQ2
 *
 * Ensure that R = scalarP * P + scalarQ * Q
 *
 * @param scalarP the computed dlog
 * @param scalarQ the computed dlog
 * @param PQ2 a 2^f-torsion basis
 * @param R a point of order dividing 2^f
 */
void ec_dlog_2(digit_t* scalarP, digit_t* scalarQ,
    const ec_basis_t* PQ2, const ec_point_t* R, const ec_curve_t* curve);

/**
 * @brief Compute the generalized dlog of R wrt the 3^e-basis PQ3
 *
 * Ensure that R = scalarP * P + scalarQ * Q
 *
 * @param scalarP the computed dlog
 * @param scalarQ the computed dlog
 * @param PQ3 a 3^e-torsion basis
 * @param R a point of order dividing 3^e
 */
void ec_dlog_3(digit_t* scalarP, digit_t* scalarQ,
    const ec_basis_t* PQ3, const ec_point_t* R, const ec_curve_t* curve);
/** @}
*/

/** @defgroup ec_isog_t Isogenies
 * @{
*/

/**
 * @brief Evaluate isogeny of even degree on list of points
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of points to evaluate the isogeny on, modified in place
 * @param length of the list points
 */
void ec_eval_even(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_point_t* points, unsigned short length);

/**
 * @brief Evaluate isogeny of even degree on list of points, assuming the point (0,0) is not in the kernel
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of points to evaluate the isogeny on, modified in place
 * @param length of the list points
 */
void ec_eval_even_nonzero(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_point_t* points, unsigned short length);

/**
 * @brief Evaluate isogeny of even degree on list of torsion bases
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of bases to evaluate the isogeny on, modified in place
 * @param length of the list bases
 */
static inline void ec_eval_even_basis(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_basis_t* points, unsigned short length) {
    ec_eval_even(image, phi, (ec_point_t*)points, sizeof(ec_basis_t) / sizeof(ec_point_t) * length);
}

/**
 * @brief Evaluate isogeny of odd degree on list of points
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of points to evaluate the isogeny on, modified in place
 * @param length of the list points
 */
//void ec_eval_odd(ec_curve_t* image, const ec_isog_odd_t* phi,
    //ec_point_t* points, unsigned short length);

//void ec_eval_three(ec_curve_t* image, const ec_isog_odd_t* phi,
        //ec_point_t* points, unsigned short length);

/**
 * @brief Evaluate isogeny of odd degree on list of torsion bases
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of bases to evaluate the isogeny on, modified in place
 * @param length of the list bases
 */
//static inline void ec_eval_odd_basis(ec_curve_t* image, const ec_isog_odd_t* phi,
    //ec_basis_t* points, unsigned short length) {
    //ec_eval_odd(image, phi, (ec_point_t*)points, sizeof(ec_basis_t) / sizeof(ec_point_t) * length);
//}

/** @}
*/

// end ec
/** @}
*/

/* Added from former curve_extras.h */

/**
 * @brief Point equals point at infinity (zero)
 *
 * @param P a point
 * @return 1 if P is the point at infinity
 */
bool ec_is_zero(ec_point_t const* P);

/**
 * @brief Sets point to point at infinity (zero)
 *
 * @param P a point
 * 
 * Sets P to point at infinity
 */
void ec_set_zero(ec_point_t *P);

/**
 * @brief Copy elliptic curve
 *
 * @param E1 an elliptic curve
 * @param E2 an elliptic curve
 * 
 * Copies E2<-E1.
 */
void copy_curve(ec_curve_t* E1, ec_curve_t const* E2);

/**
 * @brief Copy an (elliptic curve) point
 *
 * @param P a point
 * @param Q a point
 * 
 * Copies Q<-P.
 */
void copy_point(ec_point_t* P, ec_point_t const* Q);

/**
 * @brief Swap two points
 *
 * @param P a point
 * @param Q a point
 * @param option
 * 
 * If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P.
 */
void swap_points(ec_point_t* P, ec_point_t* Q, const digit_t option);

/**
 * @brief Point initialization
 *
 * @param P a point
 * 
 * P is initialized to the point at infinity.
 */
void ec_init(ec_point_t* P);

/**
* @brief Point doubling (x-only, version 1)
 *
 * @param Q a point
 * @param P a point
 * @param AC an elliptic curve
 * 
 * Sets Q<-2*P.
 */
void xDBL(ec_point_t* Q, ec_point_t const* P, ec_point_t const* AC);

/**
 * @brief Point doubling (x-only, version 2)
 *
 * @param Q a point
 * @param P a point
 * @param A24 a point representing the elliptic curve (A:C), A24=(A+2C:4C).
 * 
 * Sets Q<-2*P.
 */
void xDBLv2(ec_point_t* Q, ec_point_t const* P, ec_point_t const* A24);

/**
 * @brief Differential addition
 *
 * @param R a point
 * @param P a point
 * @param Q a point
 * @param PQ=P-Q
 * 
 * Sets R<-P+Q.
 */
void xADD(ec_point_t* R, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ);

/**
 * @brief Point doubling and addition (x-only)
 *
 * @param R a point
 * @param S a point
 * @param P a point
 * @param Q a point
 * @param PQ=P-Q
 * @param A24 a point representing the elliptic curve (A:C), A24=(A+2C:4C).
 * 
 * Sets R<-2*P and S<-2*P+Q
 */
void xDBLADD(ec_point_t* R, ec_point_t* S, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24);

/**
 * @brief Point scalar multiplication
 *
 * @param Q a point
 * @param P a point
 * @param k a multiprecision integer
 * @param curve an elliptic curve
 * 
 * Sets Q<-k*P.
 */
void xMUL(ec_point_t* Q, ec_point_t const* P, digit_t const* k, ec_curve_t const* curve);

/**
 * @brief Point linear combination (x-only)
 *
 * @param S a point
 * @param P a point
 * @param k a multiprecision integer
 * @param Q a point
 * @param l a multiprecision integer
 * @param PQ=P-Q
 * @param curve an elliptic curve
 * 
 * Sets S<-k*P+l*Q.
 */
void xDBLMUL(ec_point_t* S, ec_point_t const* P, digit_t const* k, ec_point_t const* Q, digit_t const* l, ec_point_t const* PQ, ec_curve_t const* curve);

/**
 * @brief Converting point to curve
 *
 * @param curve an elliptic curve
 * @param AC a point (A:C)
 * 
 * Sets curve coefficients to (A:C).
 */
void point_to_curve(ec_curve_t* curve, ec_point_t const * AC);

/**
 * @brief Compute (A+2C:4C) from point
 *
 * @param A24 a point
 * @param AC a point (A:C)
 * 
 * Sets A24<-(A+2C:4C).
 */
void AC_to_A24(ec_point_t* A24, ec_point_t const * AC);

/**
 * @brief Compute (A+2C:4C) from curve
 *
 * @param A24 a point
 * @param AC an elliptic curve with coefficients (A:C)
 * 
 * Sets A24<-(A+2C:4C).
 */
void curve_to_A24(ec_point_t* A24, ec_curve_t const * curve);

#define is_point_equal ec_is_equal
#define xADD ec_add

/* Jacobian model */

static void jac_init(jac_point_t* P);

bool is_jac_equal(const jac_point_t* P, const jac_point_t* Q);

static bool is_jac_xz_equal(const jac_point_t* P, const ec_point_t* Q);

static void copy_jac_point(jac_point_t* P, jac_point_t const* Q);

static void jac_neg(jac_point_t* Q, jac_point_t const* P);

void DBL(jac_point_t* Q, jac_point_t const* P, ec_curve_t const* AC);

void ADD(jac_point_t* R, jac_point_t const* P, jac_point_t const* Q, ec_curve_t const* AC);

void TPL(jac_point_t* Q, jac_point_t const* P, ec_curve_t const* AC);

void recover_y(fp2_t* y, fp2_t const* Px, ec_curve_t const* curve);

void jac_swap_points(jac_point_t* P, jac_point_t* Q, const digit_t option);

void jac_MUL(jac_point_t* Q, const jac_point_t* P, const digit_t* k, const ec_curve_t* curve);

void DBLMUL(jac_point_t* R, const jac_point_t* P, const digit_t k, const jac_point_t* Q, const digit_t l, const ec_curve_t* curve);


#endif



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////// ORIGINAL VERSION

#if 0

/** @file
 *
 * @authors Luca De Feo, Francisco RH
 *
 * @brief Elliptic curve stuff
*/

#ifndef EC_H
#define EC_H

#include <gf.h>

/** @defgroup ec Elliptic curves
 * @{
*/

/** @defgroup ec_t Data structures
 * @{
*/

/** @brief Projective point
 *
 * @typedef ec_point_t
 *
 * @struct ec_point_t
 *
 * A projective point in (X:Z) or (X:Y:Z) coordinates (tbd).
*/
typedef struct ec_point_t {
    fp2_t X;
    //fp2_t Y;
    fp2_t Z;
} ec_point_t;

/** @brief A basis of a torsion subgroup
 *
 * @typedef ec_basis_t
 *
 * @struct ec_basis_t
 *
 * A pair of points (or a triplet, tbd) forming a basis of a torsion subgroup.
*/
typedef struct ec_basis_t {
    ec_point_t P;
    ec_point_t Q;
    ec_point_t PmQ;  // or maybe not
} ec_basis_t;

/** @brief An elliptic curve
 *
 * @typedef ec_curve_t
 *
 * @struct ec_curve_t
 *
 * An elliptic curve in projective Montgomery form
*/
typedef struct ec_curve_t {
    fp2_t A;
    fp2_t C; ///< cannot be 0
} ec_curve_t;

/** @brief An isogeny of degree a power of 2
 *
 * @typedef ec_isog_even_t
 *
 * @struct ec_isog_even_t
*/
typedef struct ec_isog_even_t {
    ec_curve_t curve;      ///< The domain curve
    ec_point_t kernel;     ///< A kernel generator
    unsigned short length; ///< The length as a 2-isogeny walk
} ec_isog_even_t;


/** @brief An odd divisor of p² - 1
 *
 * @typedef ec_isog_odd_t
 *
 * @struct ec_isog_odd_t
 *
 * Given that the list of divisors of p² - 1 is known, this could be
 * represented as a fixed-length vector of integer exponents, possibly
 * distinguishing the divisors of p + 1 from those of p - 1.
*/
typedef struct ec_degree_odd_t {
    // todo (basically a ushort vector)
} ec_degree_odd_t;

/** @brief An isogeny of odd degree dividing p² - 1
 *
 * @typedef ec_isog_odd_t
 *
 * @struct ec_isog_odd_t
*/
typedef struct ec_isog_odd_t {
    ec_point_t ker_plus;    ///< A generator of E[p+1] ∩ ker(φ)
    ec_point_t ker_minus;   ///< A generator of E[p-1] ∩ ker(φ)
    ec_degree_odd_t degree; ///< The degree of the isogeny
} ec_isog_odd_t;

/** @brief Isomorphism of Montgomery curves
 *
 * @typedef ec_isom_t
 *
 * @struct ec_isom_t
 *
 * The isomorphism is given by the map maps (X:Z) ↦ ( (Nx X - Nz Z) : (D Z) )
 * TODO: fix if (X:Y:Z) coordinates.
*/
typedef struct ec_isom_t {
    fp2_t Nx;
    fp2_t Nz;
    fp2_t D;
} ec_isom_t;

// end ec_t
/** @}
*/


/** @defgroup ec_curve_t Curves and isomorphisms
 * @{
*/

/**
 * @brief j-invariant.
 *
 * @param j_inv computed j_invariant
 * @param curve input curve
 */
void ec_j_inv(fp2_t* j_inv, const ec_curve_t* curve);

/**
 * @brief Isomorphism of elliptic curve
 *
 * @param isom computed isomorphism
 * @param from domain curve
 * @param to image curve
 */
void ec_isomorphism(ec_isom_t* isom, const ec_curve_t* from, const ec_curve_t* to);

/** @}
*/
/** @defgroup ec_point_t Point operations
 * @{
*/

/**
 * @brief Point equality
 *
 * @param P a point
 * @param Q a point
 * @return 1 if equal
 */
int ec_is_equal(const ec_point_t* P, const ec_point_t* Q);

/**
 * @brief Reduce Z-coordinate of point in place
 *
 * @param P a point
 */
void ec_normalize(ec_point_t* P);

/**
 * @brief Test whether a point is on a curve
 *
 * @param curve a curve
 * @param P a point
 * @return 1 if P is on the curve
 */
int ec_is_on_curve(const ec_curve_t* curve, const ec_point_t* P);

/**
 * @brief Point negation
 *
 * @param res computed opposite of P
 * @param P a point
 */
void ec_neg(ec_point_t* res, const ec_point_t* P);

/**
 * @brief Point addition
 *
 * Needs to be adjusted if (X:Z) arithmetic.
 *
 * @param res computed sum of P and Q
 * @param P a point
 * @param Q a point
 */
void ec_add(ec_point_t* res, const ec_point_t* P, const ec_point_t* Q);

/**
 * @brief Point doubling
 *
 * @param res computed double of P
 * @param P a point
 */
void ec_dbl(ec_point_t* res, const ec_curve_t* curve, const ec_point_t* P);

/**
 * @brief Point multiplication
 *
 * @param res computed scalar * P
 * @param curve the curve
 * @param scalar an unsigned multi-precision integer
 * @param P a point
 */
void ec_mul(ec_point_t* res, const ec_curve_t* curve, const digit_t* scalar, const ec_point_t* P);

/**
 * @brief Linear combination of points of a basis
 *
 * @param res computed scalarP * P + scalarQ * Q
 * @param curve the curve
 * @param scalarP an unsigned multi-precision integer
 * @param scalarQ an unsigned multi-precision integer
 * @param PQ a torsion basis consisting of points P and Q
 */
void ec_biscalar_mul(ec_point_t* res, const ec_curve_t* curve,
    const digit_t* scalarP, const digit_t* scalarQ,
    const ec_basis_t* PQ);

/** @}
*/

/** @defgroup ec_dlog_t Discrete logs and bases
 * @{
*/

/**
 * @brief Generate a Montgomery curve and a 2^f-torsion basis
 *
 * The algorithm is deterministc
 *
 * @param curve the computed curve
 * @param PQ2 a basis of the 2^f-torsion
 * @param j_inv a j-invariant
 */
void ec_j_to_basis_2(ec_curve_t* curve, ec_basis_t* PQ2, const fp2_t* j_inv);

/**
 * @brief Complete a basis of the 2^f-torsion
 *
 * The algorithm is deterministc
 *
 * @param PQ2 a basis of the 2^f-torsion containing P as first generator
 * @param curve the curve
 * @param P a point of order 2^f
 */
void ec_complete_basis_2(ec_basis_t* PQ2, const ec_curve_t* curve, const ec_point_t* P);

/**
 * @brief Generate a 3^e-torsion basis
 *
 * The algorithm is deterministc
 *
 * @param PQ3 the computed 3^e-torsion basis
 * @param curve a curve
 */
void ec_curve_to_basis_3(ec_basis_t* PQ3, const ec_curve_t* curve);

/**
 * @brief Compute the generalized dlog of R wrt the 2^f-basis PQ2
 *
 * Ensure that R = scalarP * P + scalarQ * Q
 *
 * @param scalarP the computed dlog
 * @param scalarQ the computed dlog
 * @param PQ2 a 2^f-torsion basis
 * @param R a point of order dividing 2^f
 */
void ec_dlog_2(digit_t* scalarP, digit_t* scalarQ,
    const ec_basis_t* PQ2, const ec_point_t* R);

/**
 * @brief Compute the generalized dlog of R wrt the 3^e-basis PQ3
 *
 * Ensure that R = scalarP * P + scalarQ * Q
 *
 * @param scalarP the computed dlog
 * @param scalarQ the computed dlog
 * @param PQ3 a 3^e-torsion basis
 * @param R a point of order dividing 3^e
 */
void ec_dlog_3(digit_t* scalarP, digit_t* scalarQ,
    const ec_basis_t* PQ3, const ec_point_t* R);
/** @}
*/

/** @defgroup ec_isog_t Isogenies
 * @{
*/

/**
 * @brief Evaluate isogeny of even degree on list of points
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of points to evaluate the isogeny on, modified in place
 * @param length of the list points
 */
void ec_eval_even(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_point_t* points, unsigned short length);

/**
 * @brief Evaluate isogeny of even degree on list of torsion bases
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of bases to evaluate the isogeny on, modified in place
 * @param length of the list bases
 */
static inline void ec_eval_even_basis(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_basis_t* points, unsigned short length) {
    ec_eval_even(image, phi, points, sizeof(ec_basis_t) / sizeof(ec_point_t) * length);
}

/**
 * @brief Evaluate isogeny of odd degree on list of points
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of points to evaluate the isogeny on, modified in place
 * @param length of the list points
 */
void ec_eval_odd(ec_curve_t* image, const ec_isog_odd_t* phi,
    ec_point_t* points, unsigned short length);

/**
 * @brief Evaluate isogeny of odd degree on list of torsion bases
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of bases to evaluate the isogeny on, modified in place
 * @param length of the list bases
 */
static inline void ec_eval_odd_basis(ec_curve_t* image, const ec_isog_odd_t* phi,
    ec_basis_t* points, unsigned short length) {
    ec_eval_odd(image, phi, points, sizeof(ec_basis_t) / sizeof(ec_point_t) * length);
}

/** @}
*/

// end ec
/** @}
*/

#endif


#endif
