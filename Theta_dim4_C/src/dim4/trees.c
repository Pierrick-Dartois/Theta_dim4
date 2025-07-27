#include <matrices_mod4.h>
#include <trees.h>

void index_to_multindex_dim4(mod2_vec_4_t res,const unsigned int ind){
  for (int i = 0; i < 4; i++) {
    res[i] = (ind >> i) & 1;
  }
}

unsigned int multindex_to_index_dim4(const mod2_vec_4_t in){
  return in[0] + 2 * in[1] + 4 * in[2] + 8 * in[3];
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
  		mod4_mat_4x4_eval(&multind, D, &multind);
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
