#include <theta_isogenies_dim4.h>
#include <theta_structures_dim4.h>
#include <tree.h>
#include <change_theta_coords_dim4.h>
#include <hd.h>
#include <matrices_mod4.h>
#include <field.h>

/*** Functions for generic isogenies ***/

/** Codomain computation **/

/* Tree filling functions */

static inline unsigned int scal_prod(const unsigned int i,
                                     const unsigned int j) {
  unsigned int bit_i_k, bit_j_k, res;
  res = 0;
  for (int k = 0; k < 4; k++) {
    bit_i_k = (i >> k) & 1;
    bit_j_k = (j >> k) & 1;
    res = res ^ (bit_i_k & bit_j_k);
  }
  return res;
}

static inline void signed_sum(field_t *res, const field_t *in,
                              const unsigned int ind) {
  unsigned int s;
  field_set_zero(res);
  for (int j = 0; j < 16; j++) {
    s = scal_prod(ind, j);
    if (s) {
      field_add(res, res, &in[j]);
    } else {
      field_sub(res, res, &in[j]);
    }
  }
}

void fill_in_tree_deter(tree_t *T, const theta_point_dim4_t *theta_ker_8_sq,
                        const int *theta_index_to_kernel_index,
                        const unsigned int len_ker_8) {
  // In the generic case theta_index_to_kernel_index =
  // {-1,0,1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1},
  // since 1000 -> 1, 0100 -> 2, 0010 -> 4, 0001 -> 8,
  // -1 means the index is not in kernel
  unsigned int diff_index;
  int ker_index;

  for (int i = 0; i < T->n_children; i++) {
    diff_index = (unsigned int)(T->children[i].root - T->root);
    ker_index = theta_index_to_kernel_index[diff_index];
    signed_sum(&T->edges_num[i], theta_ker_8_sq[i].coords, T->children[i].root);
    signed_sum(&T->edges_den[i], theta_ker_8_sq[i].coords, T->root);
    fill_in_tree_deter(&T->children[i], theta_ker_8_sq,
                       theta_index_to_kernel_index, len_ker_8);
  }
}

void mul_tree_edges(tree_t *T, const field_t *parent_num,
                    const field_t *parent_den) {
  for (int i = 0; i < T->n_children; i++) {
    field_mul(&T->edges_num[i], &T->edges_num[i], parent_num);
    field_mul(&T->edges_den[i], &T->edges_den[i], parent_den);
    mul_tree_edges(&T->children[i], &T->edges_num[i], &T->edges_den[i]);
  }
}

void extract_tree_edges(field_t *num, field_t *den, unsigned int *children,
                        unsigned int *ind, const tree_t *T) {
  for (int i = 0; i < T->n_children; i++) {
    children[(*ind)] = T->children[i].root;
    num[(*ind)] = T->edges_num[i];
    den[(*ind)] = T->edges_den[i];
    (*ind)++;
    extract_tree_edges(num, den, children, ind, T->children[i]);
  }
}

/* High level isogeny codomain computation functions */

static void isogeny_compute_dim4(theta_struct_dim4_t *codomain, tree_t *T,
                                 const theta_point_dim4_t *kernel_8,
                                 const int *theta_index_to_kernel_index,
                                 const unsigned int len_ker_8,
                                 const unsigned int n_zeros) {
  theta_point_dim4_t theta_ker_8_sq[len_ker_8];
  for (int i = 0; i < len_ker_8; i++) {
    theta_dim4_sq(theta_ker_8_sq[i].coords, kernel_8[i]);
  }

  fill_in_tree_deter(T, theta_ker_8_sq, theta_index_to_kernel_index, len_ker_8);

  field_t parent_num, parent_den;
  field_set_one(&parent_num);
  field_set_one(&parent_den);
  mul_tree_edges(T, parent_num, parent_den);

  field_t num[15 - n_zeros], den[15 - n_zeros];
  unsigned int children[15 - n_zeros];
  unsigned int ind[1];
  ind[0] = 0;
  extract_tree_edges(num, den, children, ind, T);

  assert(ind[0] == 15 - n_zeros);

  field_t factor;
  field_proj_batched_inv_with_coeff(num, &factor, 15 - n_zeros);

  codomain->arith_precomp = 0;
  for (int i = 0; i < 16; i++) {
    field_set_zero(&codomain->inv_dual_null_point[i]);
  }
  field_copy(&codomain->inv_dual_null_point[T->root], &factor);
  for (int i = 0; i < 15 - n_zeros; i++) {
    field_mul(&codomain->inv_dual_null_point[children[i]], &num[i], &den[i]);
  }
}

void gen_isogeny_compute_dim4(theta_struct_dim4_t *codomain,
                              const theta_point_dim4_t *kernel_8, tree_t *T_gen,
                              const theta_index_to_kernel_index) {
  isogeny_compute_dim4(codomain, T_gen, kernel_8, theta_index_to_kernel_index,
                       4, 0);
}

/** Evaluation **/

void gen_isogeny_eval_dim4(theta_point_dim4_t *image,
                           const theta_point_dim4_t *point,
                           const theta_struct_dim4_t *codomain) {
  theta_dim4_sq(image->coords, point->coords);
  hadamard_dim4(image->coords, image->coords);
  for (int i = 0; i < 16; i++) {
    field_mul(&image->coords[i], &image->coords[i],
              &codomain->inv_dual_null_point[i]);
  }
  hadamard_dim4(image->coords, image->coords);
}

/*** Functions for gluing isogenies ***/

static void gluing_isogeny_special_simple_eval_T1_dim4(
    field_t *inv_fT1, const theta_point_dim4_t *T1,
    const theta_struct_dim4_t *codomain) {
  // When there are few zeros: U_i = 0 => U_{i^1} !=0 so that we don't need
  // T1+T2
  theta_dim4_sq(inv_fT1, T1->coords);
  hadamard_dim4(inv_fT1, inv_fT1);
  uint32_t ctrl;
  field_t inv_fT1_even[8], factor;
  for (int i = 0; i < 8; i++) {
    ctrl = field_is_zero(&codomain->inv_dual_null_point[2 * i]);
    // If U_{\chi}(0_B) = 0, chose U_{\chi\chi_1}(0_B) instead
    field_select(&factor, &codomain->inv_dual_null_point[2 * i],
                 &codomain->inv_dual_null_point[2 * i ^ 1], ctrl);
    field_select(&inv_fT1_even[i], &inv_fT1[2 * i], &inv_fT1[2 * i ^ 1], ctrl);
    field_mul(&inv_fT1_even[i], &inv_fT1_even[i], &factor);
  }

  field_proj_batched_inv(inv_fT1_even, 8);

  for (int i = 0; i < 8; i++) {
    field_copy(&inv_fT1[2 * i], &inv_fT1_even[i]);
    field_copy(&inv_fT1[2 * i ^ 1], &inv_fT1_even[i]);
  }
}

static inline unsigned int find_double_non_zero(field_t *coords) {
  unsigned int res = 0;
  for (int i = 0; i < 8; i++) {
    res = res ^ ((field_is_zero(&coords[i]) & field_is_zero(&coords[i ^ 2])) &
                 (res ^ i));
  }
  return res;
}

static void gluing_isogeny_special_eval_T1_dim4(
    field_t *inv_fT1, const theta_point_dim4_t *T1,
    const theta_point_dim4_t *T1pT2, const theta_struct_dim4_t *codomain) {
  // When there are many zeros so there exists i s.t. U_i = 0 and U_{i^1} = 0.
  // T1+T2 is needed.
  theta_dim4_sq(inv_fT1, T1->coords);
  hadamard_dim4(inv_fT1, inv_fT1);

  field_t HST1pT2[16], factor_translates, factor_non_translates;
  theta_dim4_sq(HST1pT2, T1pT2->coords);
  hadamard_dim4(HST1pT2, HST1pT2);

  // index corresponding to chi_0 such that U_{\chi_0}(0_B)\cdot
  // U_{\chi_0\chi_2}(0_B) != 0
  unsigned int ind_0 = find_double_non_zero(codomain->inv_dual_null_point);

  // factor_translates = HS_{\chi_0}(x)*1/U_{\chi_0}(0_B)
  field_mul(&factor_translates, &inv_fT1[ind_0],
            &codomain->inv_dual_null_point[ind_0]);
  // factor_non_translates = HS_{\chi_0\chi_2}(x+[2]T_2)*1/U_{\chi_0\chi_2}(0_B)
  field_mul(&factor_non_translates, &HST1pT2[ind_0 ^ 2],
            &codomain->inv_dual_null_point[ind_0 ^ 2]);

  uint32_t ctrl;
  field_t inv_fT1_even[8], factor1, factor2;
  for (int i = 0; i < 8; i++) {
    ctrl1 = field_is_zero(&codomain->inv_dual_null_point[2 * i]);
    ctrl2 = ctrl1 & field_is_zero(&codomain->inv_dual_null_point[2 * i ^ 1]);
    // If U_{\chi}(0_B) = 0, choose U_{\chi\chi_1}(0_B) instead
    field_select(&factor1, &codomain->inv_dual_null_point[2 * i],
                 &codomain->inv_dual_null_point[2 * i ^ 1], ctrl1);
    // If U_{\chi}(0_B) = 0 and U_{\chi\chi_1}(0_B)= 0, choose
    // U_{\chi\chi_2}(0_B) instead
    field_select(&factor1, &factor1, &codomain->inv_dual_null_point[2 * i ^ 2],
                 ctrl2);
    // Multiply by factor_non_translates if U_{\chi}(0_B) != 0 or
    // U_{\chi\chi_1}(0_B) != 0 and by factor_translates otherwise
    field_select(&factor2, &factor_non_translates, &factor_translates, ctrl2);
    field_mul(&factor2, &factor1, &factor2);

    field_select(&inv_fT1_even[i], &inv_fT1[2 * i], &inv_fT1[2 * i ^ 1], ctrl1);
    field_select(&inv_fT1_even[i], &inv_fT1_even[i], &inv_fT1[2 * i ^ 2],
                 ctrl2);

    field_mul(&inv_fT1_even[i], &inv_fT1_even[i], &factor2);
  }

  field_proj_batched_inv(inv_fT1_even, 8);

  for (int i = 0; i < 8; i++) {
    field_copy(&inv_fT1[2 * i], &inv_fT1_even[i]);
    field_copy(&inv_fT1[2 * i ^ 1], &inv_fT1_even[i]);
  }
}

void gluing_isogeny_surf_compute_dim4(gluing_isog_surf_dim4_t *isog, tree_t *T,
                                 const couple_theta_point_dim2_t *kernel_8,
                                 const int *theta_index_to_kernel_index,
                                 const unsigned int len_ker_8,
                                 const unsigned int n_zeros,
                                 const unsigned int simple, 
                                 const mod4_mat_4x4_t *A,
                                 const mod4_mat_4x4_t *B,
                                 const mod4_mat_4x4_t *C,
                                 const mod4_mat_4x4_t *D,
                                 const unsigned int is_neg, 
                                 const int nv_ind) {

  compute_change_theta_coords_dim4(&isog->mat_change_theta_coords, A, B, C, D, is_neg, nv_ind);

  theta_point_dim4_t kernel_8_out[len_ker_8];
  for(int i=0; i<len_ker_8; i++){
    couple_theta_point_dim2_to_theta_point_dim4(&kernel_8_out[i], &kernel_8[i],&isog->mat_change_theta_coords);
  }

  isogeny_compute_dim4(&isog->codomain, T, kernel_8_out,
                       theta_index_to_kernel_index, len_ker_8, n_zeros);
  if (simple) {
    gluing_isogeny_special_simple_eval_T1_dim4(isog->inv_fT1, &kernel_8_out[0],
                                               &isog->codomain);
  } else {
    gluing_isogeny_special_eval_T1_dim4(isog->inv_fT1, &kernel_8_out[0],
                                        &kernel_8_out[5], &isog->codomain);
  }

  copy_theta_point(&isog->T1.P1,&kernel_8[0].P1);
  copy_theta_point(&isog->T1.P2,&kernel_8[0].P2);
}
