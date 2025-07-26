#include <theta_isogenies_dim4.h>
#include <theta_structures_dim4.h>

/*** Functions for generic isogenies ***/

/** Codomain computation **/

/* Tree handling functions */

static inline void tree_init_gen(tree_t *T) {
  /* Initializing the tree:
  
  							         0000
  							          |
  				  ---------------------------
  				  |             |      |    |
  				1000           0100   0010 0001
  	        |             |      |
  	  --------------   -------   |
  		|       |    |   |     |   |
  	1100    1010 1001 0110 0101 0011
  	  |       |        |
   ------     |        |
   |    |     |        |
  1110 1101 1011      0111
   |
  1111
  
  */
  T->root = 0;
  T->n_children = 4;
  T->edges_num = malloc(4 * sizeof(field_t));
  T->edges_den = malloc(4 * sizeof(field_t));
  T->children = malloc(4 * sizeof(tree_t));  // 1000, 0100, 0010, 0001

  /* Level -1 */
  // Chlidren of 0000
  tree_t T1000;
  T1000->root = 1;
  T1000->n_children = 3;
  T1000->edges_num = malloc(3 * sizeof(field_t));
  T1000->edges_den = malloc(3 * sizeof(field_t));
  T1000->children = malloc(3 * sizeof(tree_t));  // 1100, 1010, 1001
  T->children[0] = T1000;

  tree_t T0100;
  T0100->root = 2;
  T0100->n_children = 2;
  T0100->edges_num = malloc(2 * sizeof(field_t));
  T0100->edges_den = malloc(2 * sizeof(field_t));
  T0100->children = malloc(2 * sizeof(tree_t));  // 0110, 0101
  T->children[1] = T0100;

  tree_t T0010;
  T0010->root = 4;
  T0010->n_children = 1;
  T0010->edges_num = malloc(1 * sizeof(field_t));
  T0010->edges_den = malloc(1 * sizeof(field_t));
  T0010->children = malloc(1 * sizeof(tree_t));  // 0011
  T->children[2] = T0010;

  tree_t T0001;
  T0001->root = 8;
  T0001->n_children = 0;
  // Terminal node
  T->children[3] = T0001;

  /* Level -2 */
  // Chlidren of 1000
  tree_t T1100;
  T1100->root = 3;
  T1100->n_children = 2;
  T1100->edges_num = malloc(2 * sizeof(field_t));
  T1100->edges_den = malloc(2 * sizeof(field_t));
  T1100->children = malloc(2 * sizeof(tree_t));  // 1110, 1101
  T1000->children[0] = T1100;

  tree_t T1010;
  T1010->root = 5;
  T1010->n_children = 1;
  T1010->edges_num = malloc(1 * sizeof(field_t));
  T1010->edges_den = malloc(1 * sizeof(field_t));
  T1010->children = malloc(1 * sizeof(tree_t));  // 1011
  T1000->children[1] = T1010;

  tree_t T1001;
  T1001->root = 9;
  T1001->n_children = 0;
  // Terminal node
  T1000->children[2] = T1001;

  // Children of 0100
  tree_t T0110;
  T0110->root = 6;
  T0110->n_children = 1;
  T0110->edges_num = malloc(1 * sizeof(field_t));
  T0110->edges_den = malloc(1 * sizeof(field_t));
  T0110->children = malloc(1 * sizeof(tree_t));  // 0111
  T0100->children[0] = T0110;

  tree_t T0101;
  T0101->root = 10;
  T0101->n_children = 0;
  // Terminal node
  T0100->children[1] = T0101;

  // Children of 0010
  tree_t T0011;
  T0011->root = 12;
  T0011->n_children = 0;
  // Terminal node
  T0010->children[0] = T0011;

  /* Level -3 */
  // Chlidren of 1100
  tree_t T1110;
  T1110->root = 7;
  T1110->n_children = 1;
  T1110->edges_num = malloc(1 * sizeof(field_t));
  T1110->edges_den = malloc(1 * sizeof(field_t));
  T1110->children = malloc(1 * sizeof(tree_t));  // 1111
  T1100->children[0] = T1110;

  tree_t T1101;
  T1101->root = 11;
  T1101->n_children = 0;
  // Terminal node
  T1100->children[1] = T1101;

  // Chlidren of 1010
  tree_t T1011;
  T1011->root = 13;
  T1011->n_children = 0;
  // Terminal node
  T1010->children[0] = T1011;

  // Chlidren of 01100
  tree_t T0111;
  T0111->root = 14;
  T0111->n_children = 0;
  // Terminal node
  T0110->children[0] = T0111;

  /* Level -4 */
  // Children of 1110
  tree_t T1111;
  T1111->root = 15;
  T1111->n_children = 0;
  // Terminal node
  T1110->children[0] = T1111;
}

void tree_finalize(tree_t *T) {
  if (T->n_children > 0) {
    free(T->edges_num);
    free(T->edges_den);
    for (int i = 0; i < T->n_children; i++) {
      tree_finalize(T->children[i]);
    }
    free(T->children);
  }
}

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

void extract_tree_edges(field_t *num, field_t *den, unsigned int *children, unsigned int *ind, const tree_t *T){
	for(int i=0; i<T->n_children; i++){
		children[(*ind)] = T->children[i].root;
		num[(*ind)] = T->edges_num[i];
		den[(*ind)] = T->edges_den[i];
		(*ind)++;
		extract_tree_edges(num, den, children, ind, T->children[i]);
	}
}

/* High level isogeny codomain computation functions */

void isogeny_compute_dim4(theta_struct_dim4_t *codomain, tree_t *T,
                          const theta_point_dim4_t *kernel_8,
                          const int *theta_index_to_kernel_index,
                          const unsigned int len_ker_8) {
	
	theta_point_dim4_t theta_ker_8_sq[len_ker_8];
  for (int i = 0; i < len_ker_8; i++) {
  	theta_dim4_sq(theta_ker_8_sq[i].coords,kernel_8[i]);
  }

  fill_in_tree_deter(T, theta_ker_8_sq, theta_index_to_kernel_index,len_ker_8);
  
  field_t parent_num, parent_den;
  field_set_one(&parent_num);
  field_set_one(&parent_den);
  mul_tree_edges(T, parent_num, parent_den);

  field_t num[15], den[15];
  unsigned int children[15];
  unsigned int ind[1];
  ind[0] = 0;
  extract_tree_edges(num, den, children, ind, T);

  //assert(ind[0] == 15);

  field_t factor;
  fp_proj_batched_inv_with_coeff(num,&factor,15);

  codomain->arith_precomp = 0;
  field_copy(&codomain->inv_dual_null_point[T->root],&factor);
  for(int i=0; i<15; i++){
  	field_mul(&codomain->inv_dual_null_point[children[i]],&num[i],&den[i]);
  }
}

void gen_isogeny_compute_dim4(theta_struct_dim4_t *codomain, const theta_point_dim4_t *kernel_8){
	tree_t T;
	tree_init_gen(&T);
	int theta_index_to_kernel_index[16] = {-1,0,1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1};

	isogeny_compute_dim4(codomain, &T, kernel_8, theta_index_to_kernel_index,4);
	tree_finalize(&T);
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



