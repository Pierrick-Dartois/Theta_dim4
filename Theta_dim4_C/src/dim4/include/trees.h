#ifndef TREES_H
#define TREES_H

#include <matrices_mod4.h>
#include <theta_structures.h>
#include <field.h>

typedef struct tree {
  unsigned int len_support;
  int *support;
  int[16][16] adj_mat;
  field_t *edges_parents;
  field_t *edges_children;

  int n_children;
  tree_t *children;
  field_t *edges_num;
  field_t *edges_den;
} tree_t;

/* Deprecated

typedef struct tree {
  int root;
  int n_children;
  tree_t *children;
  field_t *edges_num;
  field_t *edges_den;
} tree_t;

void tree_init_gen(tree_t *T);
void find_zeros_surf(int *zeros, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D);
void find_zeros_EC(int *zeros, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D);
unsigned int fill_in_adj_matrix(unsigned int **mat,  const int *zeros);
void build_tree_from_adj_matrix(tree_t *T, const unsigned int **mat, const int root);
void tree_finalize(tree_t *T);
*/

#endif