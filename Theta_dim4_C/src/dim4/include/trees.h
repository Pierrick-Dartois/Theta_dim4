#ifndef TREES_H
#define TREES_H

#include <matrices_mod4.h>
#include <theta_structures.h>

typedef struct tree {
  int root;
  int n_children;
  tree_t *children;
  field_t *edges_num;
  field_t *edges_den;
} tree_t;

void tree_init_gen(tree_t *T);
void index_to_multindex_dim4(mod2_vec_4_t res,const unsigned int ind);
unsigned int multindex_to_index_dim4(const mod2_vec_4_t in);
void find_zeros(int *zeros, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D);
unsigned int fill_in_adj_matrix(unsigned int **mat,  const int *zeros);
void build_tree_from_adj_matrix(tree_t *T, const unsigned int **mat, const int root);
void tree_finalize(tree_t *T);

#endif