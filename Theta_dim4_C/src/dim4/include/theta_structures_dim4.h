#ifndef THETA_STRUCTURES_DIM4_H
#define THETA_STRUCTURES_DIM4_H

// Relevant includes: TODO
#include <field.h>
#include <theta_structure.h>//from dim2

typedef struct theta_point_dim4 {
  field_t coords[16];
} theta_point_dim4_t;

typedef struct theta_struct_dim4 {
  field_t inv_dual_null_point[16];  // Stores the inverse dual theta null point
                                    // by default (for isogeny evaluations).
                                    // This saves a batch inversion for every
                                    // isogeny computation
  int arith_precomp;
  theta_point_dim4_t null_point;
  filed_t inv_codomain_dual_null_point_sq[16];
} theta_struct_dim4_t;

typedef struct couple_theta_point_dim2 {
  theta_point_t P1;
  theta_point_t P2;
} couple_theta_point_dim2_t;

typedef struct couple_theta_struct_dim2 {
  theta_structure_t A1;
  theta_structure_t A2;
} couple_theta_struct_dim2_t;

void hadamard_dim4(field_t *out, const field_t *in);
void theta_dim4_sq(field_t *out, const field_t *in);
void theta_struct_arith_precomp(theta_struct_dim4_t *theta_struct);
void theta_double(theta_point_dim4_t *out, const theta_point_dim4_t *in,
                  theta_struct_dim4_t *theta_struct);
void theta_double_iter(theta_point *out, const theta_point_dim4_t *in,
                       theta_struct_dim4_t *theta_struct, int n);
void copy_theta_point_dim4(theta_point_dim4_t *out, 
  const theta_point_dim4_t *in);

#endif