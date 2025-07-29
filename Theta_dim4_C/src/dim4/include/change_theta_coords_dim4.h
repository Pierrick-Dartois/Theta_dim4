#ifndef CHANGE_THETA_COORDS_DIM4_H
#define CHANGE_THETA_COORDS_DIM4_H

#include <field.h>
#include <matrices_mod4.h>
#include <theta_structure.h>  // Imported from dim2
#include <theta_structures_dim4.h>

typedef field_t field_mat_16x16_t[16][16];

void compute_change_theta_coords_dim4(field_mat_16x16_t *mat,
                                      const mod4_mat_4x4_t *A,
                                      const mod4_mat_4x4_t *B,
                                      const mod4_mat_4x4_t *C,
                                      const mod4_mat_4x4_t *D, 
                                      const unsigned int is_neg, 
                                      const int nv_ind);

void couple_theta_point_dim2_to_theta_point_dim4(
    theta_point_dim4_t *res, const couple_theta_point_dim2_t *in,
    const field_mat_16x16_t *mat_change_theta_coords);

#endif