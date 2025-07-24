#include <theta_structures_dim4.h>

typedef struct tree {
  int root;
  int n_children;
  tree_t *children;
  field_t *edges_num;
  field_t *edges_den;
} tree_t;

void gen_isogeny_compute_dim4(theta_struct_dim4_t *codomain, const theta_point_dim4_t *kernel_8);
void gen_isogeny_eval_dim4(theta_point_dim4_t *image, const theta_point_dim4_t *point, const theta_struct_dim4_t *codomain);