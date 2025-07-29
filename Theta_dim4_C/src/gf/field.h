#include <fp.h>
#include <fp2.h>

#ifdef FP_ONLY

#define field_t fp_t
#define field_set_small fp_set_small
#define field_mul_small fp_mul_small
#define field_set_zero fp_set_zero
#define field_set_one fp_set_one
#define field_is_equal fp_is_equal
#define field_is_zero fp_is_zero
#define field_is_one fp_is_one
#define field_copy fp_copy
#define field_encode fp_encode
#define field_decode fp_decode
#define field_select fp_select
#define field_cswap fp_cswap
#define field_add fp_add
#define field_sub fp_sub
#define field_neg fp_neg
#define field_mul fp_mul
#define field_sqr fp_sqr
#define field_inv fp_inv
#define field_is_square fp_is_square
#define field_sqrt fp_sqrt
#define field_half fp_half
// void fp_exp3div4(fp_t *out, const fp_t *a);
// void fp_div3(fp_t *out, const fp_t *a);
#define field_batched_inv fp_batched_inv
#define field_proj_batched_inv fp_proj_batched_inv
#define field_proj_batched_inv_with_coeff fp_proj_batched_inv_with_coeff
#define field_print fp_print

#else

#define field_t fp2_t
#define field_set_small fp2_set_small
#define field_mul_small fp2_mul_small
#define field_set_zero fp2_set_zero
#define field_set_one fp2_set_one
#define field_is_equal fp2_is_equal
#define field_is_zero fp2_is_zero
#define field_is_one fp2_is_one
#define field_copy fp2_copy
#define field_encode fp2_encode
#define field_decode fp2_decode
#define field_select fp2_select
#define field_cswap fp2_cswap
#define field_add fp2_add
// void fp2_add_one(fp2_t *x, const fp2_t *y);
#define field_sub fp2_sub
#define field_neg fp2_neg
#define field_mul fp2_mul
#define field_sqr fp2_sqr
#define field_inv fp2_inv
#define field_is_square fp2_is_square
#define field_sqrt fp2_sqrt
// uint32_t fp2_sqrt_verify(fp2_t *a);
#define field_half fp2_half
#define field_batched_inv fp2_batched_inv
#define field_proj_batched_inv fp2_proj_batched_inv
#define field_proj_batched_inv_with_coeff fp2_proj_batched_inv_with_coeff
// void fp2_pow_vartime(fp2_t *out, const fp2_t *x, const digit_t *exp, const
// int size);
#define field_print fp2_print

#endif