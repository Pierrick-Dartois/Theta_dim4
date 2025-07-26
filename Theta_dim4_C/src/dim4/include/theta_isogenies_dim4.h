#ifndef THETA_ISOGENIES_DIM4_H
#define THETA_ISOGENIES_DIM4_H

#include <theta_structures_dim4.h>
#include <trees.h>


void gen_isogeny_compute_dim4(theta_struct_dim4_t *codomain, const theta_point_dim4_t *kernel_8);
void gen_isogeny_eval_dim4(theta_point_dim4_t *image, const theta_point_dim4_t *point, const theta_struct_dim4_t *codomain);

#endif