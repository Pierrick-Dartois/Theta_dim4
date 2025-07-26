#ifndef TREES_H
#define TREES_H

#include <matrices_mod4.h>

typedef struct tree {
  int root;
  int n_children;
  tree_t *children;
  field_t *edges_num;
  field_t *edges_den;
} tree_t;

void find_zeros(unsigned int *zero_index, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D);
void fill_in_adj_matrix(unsigned int **mat, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D);

#endif