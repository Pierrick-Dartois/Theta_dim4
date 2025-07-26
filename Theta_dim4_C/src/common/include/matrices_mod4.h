#ifndef MATRICES_MOD4_H
#define MATRICES_MOD4_H

extern const unsigned int INV_MOD4[4]={0,1,0,3};

typedef unsigned int mod2_vec_2_t[2];

typedef unsigned int mod2_vec_4_t[4];

typedef unsigned int mod4_mat_2x2_t[2][2];

typedef unsigned int mod4_mat_4x4_t[4][4];

typedef unsigned int mod4_mat_8x8_t[8][8];

void mod4_mat_2x2_eval(mod2_vec_2_t *res, const mod4_mat_2x2_t *mat, const mod2_vec_2_t *vec);
void mod4_mat_2x2_copy(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *mat);
void mod4_mat_2x2_mul(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *A, const mod4_mat_2x2_t *B);
void mod4_mat_2x2_transpose(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *mat);
uint32_t mod4_mat_2x2_inv(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *mat);

void mod4_mat_4x4_eval(mod2_vec_4_t *res, const mod4_mat_4x4_t *mat, const mod2_vec_4_t *vec);
void mod4_mat_4x4_copy(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *mat);
void mod4_mat_4x4_mul(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *A, const mod4_mat_4x4_t *B);
void mod4_mat_4x4_transpose(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *mat);
uint32_t mod4_mat_4x4_inv(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *mat);

void mod4_scal_prod_2(unsigned int res, const mod2_vec_2_t *a, const mod2_vec_2_t *b);
void mod2_scal_prod_2(unsigned int res, const mod2_vec_2_t *a, const mod2_vec_2_t *b);

void mod4_scal_prod_4(unsigned int res, const mod2_vec_4_t *a, const mod2_vec_4_t *b);
void mod2_scal_prod_4(unsigned int res, const mod2_vec_4_t *a, const mod2_vec_4_t *b);

#endif

