#include <matrices_mod4.h>
#include <trees.h>

void tree_init_gen(tree_t *T) {
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

void index_to_multindex_dim4(mod2_vec_4_t res,const unsigned int ind){
  for (int i = 0; i < 4; i++) {
    res[i] = (ind >> i) & 1;
  }
}

unsigned int multindex_to_index_dim4(const mod2_vec_4_t in){
  return (in[0]&1) + 2 * (in[1]&1) + 4 * (in[2]&1) + 8 * (in[3]&1);
}

void find_zeros(int *zeros, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D) {
  mod4_mat_4x4_t tDC;

  mod4_mat_4x4_transpose(&tDC, D);
  mod4_mat_4x4_mul(&tDC, &tD, C);

  int c = (tDC[0][2] & 1) + 2 * (tDC[0][3] & 1) + 4 * (tDC[1][2] & 1) +
          8 * (tDC[1][3] & 1);

  // Finding the zeros of the normalised form [[I4, tDC],[0, I4]].
  switch (c) {
    case 1:
      zeros = {5, 7, 13, 15, -1, -1};
      break;
    case 2:
      zeros = {9, 11, 13, 15, -1, -1};
      break;
    case 3:
      zeros = {5, 7, 9, 11, -1, -1};
      break;
    case 4:
      zeros = {6, 7, 14, 15, -1, -1};
      break;
    case 5:
      zeros = {5, 6, 13, 14, -1, -1};
      break;
    case 6:
      zeros = {6, 7, 9, 11, 13, 14};
      break;
    case 7:
      zeros = {5, 6, 9, 11, 14, 15};
      break;
    case 8:
      zeros = {10, 11, 14, 15, -1, -1};
      break;
    case 9:
      zeros = {5, 7, 10, 11, 13, 14};
      break;
    case 10:
      zeros = {9, 10, 13, 14, -1, -1};
      break;
    case 11:
      zeros = {5, 7, 9, 10, 14, 15};
      break;
    case 12:
      zeros = {6, 7, 10, 11, -1, -1}; 
      break; 
    case 13:
      zeros = {5, 6, 10, 11, 13, 15};
      break;
    case 14:
      zeros = {6, 7, 9, 10, 13, 15};
      break;
    case 15:
      zeros = {5, 6, 9, 10, -1, -1};
      break;
    default:
      zeros = {-1, -1, -1, -1, -1, -1};
  }

  mod2_vec_4_t multind;

  // Translating the zero values (depending on D).
  for (int i=0; i<6; i++){
  	if (zeros[i]!=-1){
  		index_to_multindex_dim4(multind,zeros[i]);
  		mod2_mat_4x4_eval(&multind, D, &multind);
  		zeros[i] = multindex_to_index_dim4(multind);
  	}
  }
}

static inline uint32_t is_in(const int ind, const int *indices, const unsigned int len){
	uint32_t res;
	for(int i=0; i<len; i++){
		res = res || (ind == indices[i]);
	}
	return res;
}

static void fill_in_adj_matrix_recurse(unsigned int **mat, const int *zeros, const int ind, const int root, const int additional_point){
	// additional_point: boolean encoding if we can use T1+T2 corresponding to a translation by 1100 or not.
	if (!is_in(ind,zeros,6)) {
	int ind_trans, s;
	for(int i=0; i<4; i++){
		ind_trans = ind^(1 << i);
		s = 0;
		for(int j=0; j<16; j++){
			s = s+mat[j][ind_trans];
		}
		if ((!is_in(ind_trans,zeros,6))&&(s==0)&&(ind_trans!=root)){
			mat[ind][ind_trans] = 1;
			fill_in_adj_matrix_recurse(mat,zeros,ind_trans,root,additional_point);
		}
	}
	if(additional_point){
		ind_trans = ind^3;
		s = 0;
		for(int j=0; j<16; j++){
			s = s+mat[j][ind_trans];
		}
		if ((!is_in(ind_trans,zeros,6))&&(s==0)&&(ind_trans!=root)){
			mat[ind][ind_trans] = 1;
			fill_in_adj_matrix_recurse(mat,zeros,ind_trans,root,additional_point);
		}
	}	
	}
}


unsigned int fill_in_adj_matrix(unsigned int **mat, const int *zeros){
	for(int i=0; i<16; i++){
		for(int j=0; j<16; j++){
			mat[i][j]=0;
		}
	}

	unsigned int root = 0;
	while(is_in(root,zeros,6)){
		root += 1;
	}

	fill_in_adj_matrix_recurse(mat,zeros,root,root,0);

	// Fill in the matrix with an additional point if it is incomplete
	for(int i=0; i<16; i++){
		fill_in_adj_matrix_recurse(mat,zeros,i,root,1);
	}

	return root;
}

void build_tree_from_adj_matrix(tree_t *T, const unsigned int **mat, const int root){
  T->root = root;
  T->n_children = 0;

  for(int j=0; j<16; j++){
    T->n_children + = mat[root][j];
  }

  T->children = malloc((T->n_children)*sizeof(tree_t));
  T->edges_num = malloc((T->n_children)*sizeof(field_t));
  T->edges_den = malloc((T->n_children)*sizeof(field_t));

  int i_child = 0;
  for(int j=0; j<16; j++){
    if (mat[root][j]){
      build_tree_from_adj_matrix(T->children[i_child], mat, j);
      i_child += 1;
    }
  }
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