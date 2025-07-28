#ifndef THETA_ISOGENIES_DIM4_H
#define THETA_ISOGENIES_DIM4_H

#include <theta_structures_dim4.h>
#include <trees.h>

typedef struct gluing_isog_dim4 {
  theta_struct_dim4_t codomain;
  field_t inv_fT1[16];
} gluing_isog_dim4_t;

void gen_isogeny_compute_dim4(theta_struct_dim4_t *codomain,
                              const theta_point_dim4_t *kernel_8, tree_t *T_gen,
                              const theta_index_to_kernel_index);
void gen_isogeny_eval_dim4(theta_point_dim4_t *image,
                           const theta_point_dim4_t *point,
                           const theta_struct_dim4_t *codomain);

#endif