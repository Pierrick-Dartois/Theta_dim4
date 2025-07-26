#include <matrices_mod4.h>

void mod4_mat_2x2_eval(mod2_vec_2_t *res, const mod4_mat_2x2_t *mat,
                       const mod2_vec_2_t *vec) {
  for (int i = 0; i < 2; i++) {
    (*res)[i] = 0;
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      (*res)[i] = (*res)[i] + (*mat)[i][j] * (*vec)[j];
    }
  }
  for (int i = 0; i < 2; i++) {
    (*res)[i] = (*res)[i] & 1;
  }
}

void mod4_mat_2x2_copy(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *mat) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      (*res)[i][j] = (*mat)[i][j];
    }
  }
}

void mod4_mat_2x2_mul(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *A,
                      const mod4_mat_2x2_t *B) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      (*res)[i][j] = 0;
      for (int k = 0; k < 2; k++) {
        (*res)[i][j] = (*res)[i][j] + (*A)[i][k] * (*B)[k][j];
      }
      (*res)[i][j] = (*res)[i][j] & 0x3;
    }
  }
}

void mod4_mat_2x2_transpose(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *mat) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      (*res)[i][j] = (*mat)[j][i];
    }
  }
}

uint32_t mod4_mat_2x2_inv(mod4_mat_2x2_t *res, const mod4_mat_2x2_t *mat) {
  (*res)[0][0] = (*mat)[1][1];
  (*res)[0][1] = (-(*mat)[0][1]) & 0x3;
  (*res)[1][0] = (-(*mat)[1][0]) & 0x3;
  (*res)[1][1] = (*mat)[0][0];

  unsigned int det =
      ((*mat)[0][0] * (*mat)[1][1] - (*mat)[0][1] * (*mat)[1][0]) & 0x3;
  det = INV_MOD4[det];
  unsigned int ret = (det != 0);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      (*res)[i][j] = (*res)[i][j] & 0x3;
    }
  }
  return ret;
}

void mod4_mat_4x4_eval(mod2_vec_4_t *res, const mod4_mat_4x4_t *mat,
                       const mod2_vec_4_t *vec) {
  for (int i = 0; i < 4; i++) {
    (*res)[i] = 0;
  }
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*res)[i] = (*res)[i] + (*mat)[i][j] * (*vec)[j];
    }
  }
  for (int i = 0; i < 4; i++) {
    (*res)[i] = (*res)[i] & 1;
  }
}

void mod4_mat_4x4_copy(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *mat) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*res)[i][j] = (*mat)[i][j];
    }
  }
}

void mod4_mat_4x4_mul(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *A,
                      const mod4_mat_4x4_t *B) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*res)[i][j] = 0;
      for (int k = 0; k < 4; k++) {
        (*res)[i][j] = (*res)[i][j] + (*A)[i][k] * (*B)[k][j];
      }
      (*res)[i][j] = (*res)[i][j] & 0x3;
    }
  }
}

void mod4_mat_4x4_transpose(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *mat) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*res)[i][j] = (*mat)[j][i];
    }
  }
}

uint32_t mod4_mat_4x4_inv(mod4_mat_4x4_t *res, const mod4_mat_4x4_t *mat) {
  (*res)[0][0] = (*mat)[1][1] * (*mat)[2][2] * (*mat)[3][3] -
                 (*mat)[1][1] * (*mat)[2][3] * (*mat)[3][2] -
                 (*mat)[2][1] * (*mat)[1][2] * (*mat)[3][3] +
                 (*mat)[2][1] * (*mat)[1][3] * (*mat)[3][2] +
                 (*mat)[3][1] * (*mat)[1][2] * (*mat)[2][3] -
                 (*mat)[3][1] * (*mat)[1][3] * (*mat)[2][2];

  (*res)[1][0] = -(*mat)[1][0] * (*mat)[2][2] * (*mat)[3][3] +
                 (*mat)[1][0] * (*mat)[2][3] * (*mat)[3][2] +
                 (*mat)[2][0] * (*mat)[1][2] * (*mat)[3][3] -
                 (*mat)[2][0] * (*mat)[1][3] * (*mat)[3][2] -
                 (*mat)[3][0] * (*mat)[1][2] * (*mat)[2][3] +
                 (*mat)[3][0] * (*mat)[1][3] * (*mat)[2][2];

  (*res)[2][0] = (*mat)[1][0] * (*mat)[2][1] * (*mat)[3][3] -
                 (*mat)[1][0] * (*mat)[2][3] * (*mat)[3][1] -
                 (*mat)[2][0] * (*mat)[1][1] * (*mat)[3][3] +
                 (*mat)[2][0] * (*mat)[1][3] * (*mat)[3][1] +
                 (*mat)[3][0] * (*mat)[1][1] * (*mat)[2][3] -
                 (*mat)[3][0] * (*mat)[1][3] * (*mat)[2][1];

  (*res)[3][0] = -(*mat)[1][0] * (*mat)[2][1] * (*mat)[3][2] +
                 (*mat)[1][0] * (*mat)[2][2] * (*mat)[3][1] +
                 (*mat)[2][0] * (*mat)[1][1] * (*mat)[3][2] -
                 (*mat)[2][0] * (*mat)[1][2] * (*mat)[3][1] -
                 (*mat)[3][0] * (*mat)[1][1] * (*mat)[2][2] +
                 (*mat)[3][0] * (*mat)[1][2] * (*mat)[2][1];

  (*res)[0][1] = -(*mat)[0][1] * (*mat)[2][2] * (*mat)[3][3] +
                 (*mat)[0][1] * (*mat)[2][3] * (*mat)[3][2] +
                 (*mat)[2][1] * (*mat)[0][2] * (*mat)[3][3] -
                 (*mat)[2][1] * (*mat)[0][3] * (*mat)[3][2] -
                 (*mat)[3][1] * (*mat)[0][2] * (*mat)[2][3] +
                 (*mat)[3][1] * (*mat)[0][3] * (*mat)[2][2];

  (*res)[1][1] = (*mat)[0][0] * (*mat)[2][2] * (*mat)[3][3] -
                 (*mat)[0][0] * (*mat)[2][3] * (*mat)[3][2] -
                 (*mat)[2][0] * (*mat)[0][2] * (*mat)[3][3] +
                 (*mat)[2][0] * (*mat)[0][3] * (*mat)[3][2] +
                 (*mat)[3][0] * (*mat)[0][2] * (*mat)[2][3] -
                 (*mat)[3][0] * (*mat)[0][3] * (*mat)[2][2];

  (*res)[2][1] = -(*mat)[0][0] * (*mat)[2][1] * (*mat)[3][3] +
                 (*mat)[0][0] * (*mat)[2][3] * (*mat)[3][1] +
                 (*mat)[2][0] * (*mat)[0][1] * (*mat)[3][3] -
                 (*mat)[2][0] * (*mat)[0][3] * (*mat)[3][1] -
                 (*mat)[3][0] * (*mat)[0][1] * (*mat)[2][3] +
                 (*mat)[3][0] * (*mat)[0][3] * (*mat)[2][1];

  (*res)[3][1] = (*mat)[0][0] * (*mat)[2][1] * (*mat)[3][2] -
                 (*mat)[0][0] * (*mat)[2][2] * (*mat)[3][1] -
                 (*mat)[2][0] * (*mat)[0][1] * (*mat)[3][2] +
                 (*mat)[2][0] * (*mat)[0][2] * (*mat)[3][1] +
                 (*mat)[3][0] * (*mat)[0][1] * (*mat)[2][2] -
                 (*mat)[3][0] * (*mat)[0][2] * (*mat)[2][1];

  (*res)[0][2] = (*mat)[0][1] * (*mat)[1][2] * (*mat)[3][3] -
                 (*mat)[0][1] * (*mat)[1][3] * (*mat)[3][2] -
                 (*mat)[1][1] * (*mat)[0][2] * (*mat)[3][3] +
                 (*mat)[1][1] * (*mat)[0][3] * (*mat)[3][2] +
                 (*mat)[3][1] * (*mat)[0][2] * (*mat)[1][3] -
                 (*mat)[3][1] * (*mat)[0][3] * (*mat)[1][2];

  (*res)[1][2] = -(*mat)[0][0] * (*mat)[1][2] * (*mat)[3][3] +
                 (*mat)[0][0] * (*mat)[1][3] * (*mat)[3][2] +
                 (*mat)[1][0] * (*mat)[0][2] * (*mat)[3][3] -
                 (*mat)[1][0] * (*mat)[0][3] * (*mat)[3][2] -
                 (*mat)[3][0] * (*mat)[0][2] * (*mat)[1][3] +
                 (*mat)[3][0] * (*mat)[0][3] * (*mat)[1][2];

  (*res)[2][2] = (*mat)[0][0] * (*mat)[1][1] * (*mat)[3][3] -
                 (*mat)[0][0] * (*mat)[1][3] * (*mat)[3][1] -
                 (*mat)[1][0] * (*mat)[0][1] * (*mat)[3][3] +
                 (*mat)[1][0] * (*mat)[0][3] * (*mat)[3][1] +
                 (*mat)[3][0] * (*mat)[0][1] * (*mat)[1][3] -
                 (*mat)[3][0] * (*mat)[0][3] * (*mat)[1][1];

  (*res)[3][2] = -(*mat)[0][0] * (*mat)[1][1] * (*mat)[3][2] +
                 (*mat)[0][0] * (*mat)[1][2] * (*mat)[3][1] +
                 (*mat)[1][0] * (*mat)[0][1] * (*mat)[3][2] -
                 (*mat)[1][0] * (*mat)[0][2] * (*mat)[3][1] -
                 (*mat)[3][0] * (*mat)[0][1] * (*mat)[1][2] +
                 (*mat)[3][0] * (*mat)[0][2] * (*mat)[1][1];

  (*res)[0][3] = -(*mat)[0][1] * (*mat)[1][2] * (*mat)[2][3] +
                 (*mat)[0][1] * (*mat)[1][3] * (*mat)[2][2] +
                 (*mat)[1][1] * (*mat)[0][2] * (*mat)[2][3] -
                 (*mat)[1][1] * (*mat)[0][3] * (*mat)[2][2] -
                 (*mat)[2][1] * (*mat)[0][2] * (*mat)[1][3] +
                 (*mat)[2][1] * (*mat)[0][3] * (*mat)[1][2];

  (*res)[1][3] = (*mat)[0][0] * (*mat)[1][2] * (*mat)[2][3] -
                 (*mat)[0][0] * (*mat)[1][3] * (*mat)[2][2] -
                 (*mat)[1][0] * (*mat)[0][2] * (*mat)[2][3] +
                 (*mat)[1][0] * (*mat)[0][3] * (*mat)[2][2] +
                 (*mat)[2][0] * (*mat)[0][2] * (*mat)[1][3] -
                 (*mat)[2][0] * (*mat)[0][3] * (*mat)[1][2];

  (*res)[2][3] = -(*mat)[0][0] * (*mat)[1][1] * (*mat)[2][3] +
                 (*mat)[0][0] * (*mat)[1][3] * (*mat)[2][1] +
                 (*mat)[1][0] * (*mat)[0][1] * (*mat)[2][3] -
                 (*mat)[1][0] * (*mat)[0][3] * (*mat)[2][1] -
                 (*mat)[2][0] * (*mat)[0][1] * (*mat)[1][3] +
                 (*mat)[2][0] * (*mat)[0][3] * (*mat)[1][1];

  (*res)[3][3] = (*mat)[0][0] * (*mat)[1][1] * (*mat)[2][2] -
                 (*mat)[0][0] * (*mat)[1][2] * (*mat)[2][1] -
                 (*mat)[1][0] * (*mat)[0][1] * (*mat)[2][2] +
                 (*mat)[1][0] * (*mat)[0][2] * (*mat)[2][1] +
                 (*mat)[2][0] * (*mat)[0][1] * (*mat)[1][2] -
                 (*mat)[2][0] * (*mat)[0][2] * (*mat)[1][1];

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*res)[i][j] = (*res)[i][j] & 0x3;  // Reduction mod 4
    }
  }

  unsigned int det =
      ((*mat)[0][0] * (*res)[0][0] + (*mat)[0][1] * (*res)[1][0] +
       (*mat)[0][2] * (*res)[2][0] + (*mat)[0][3] * (*res)[3][0]) &
      0x3;

  det = INV_MOD4[det];

  uint32_t ret = (det != 0);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*res)[i][j] = (*res)[i][j] * det;
    }
  }

  return ret;
}

