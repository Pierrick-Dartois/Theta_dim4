#include <change_theta_coords_dim4.h>
#include <matrices_mod4.h>
#include <theta_structure.h>
#include <theta_structures_dim4.h>
#include <fp.h>
#include <field.h>

static inline unsigned int choose_non_vanishing_index(
    field_t *row, const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D,
    const unsigned int is_neg) {
  // Should not be used for constant time applications

#ifdef FP_ONLY
  field_t one, mone;
  field_set_one(&one);
  field_neg(&mone, &one);
#else
  field_t one, mone, onexi, monexi;
  field_set_one(&one);
  field_neg(&mone, &one);
  fp_set_zero(&onexi.re);
  fp_set_one(&onexi.im);
  fp_set_zero(&monexi.re);
  fp_neg(&monexi.im, &onexi.im);
  field_cswap(&onexi, &monexi, is_neg);
#endif

  mod2_vec_4_t I, J, CJ, DJ;
  unsigned int scal, cj;
  for (int i = 0; i < 16; i++) {
    index_to_multindex_dim4(I, i);

    for (int j = 0; j < 16; j++) {
      field_set_zero(&row[j]);
    }

    for (int j = 0; j < 16; j++) {
      index_to_multindex_dim4(J, j);
      mod4_mat_4x4_eval(&CJ, C, &J);
      mod4_mat_4x4_eval(&DJ, D, &J);
      scal = (-mod4_scal_prod_4(&CJ, &DJ) - 2 * mod4_scal_prod_4(&I, &DJ)) & 3;

      cj = multindex_to_index_dim4(CJ);

#ifdef FP_ONLY

      if (scal == 0) {
        field_add(&row[cj ^ i], &row[cj ^ i], &one);
      } else {
        field_add(&row[cj ^ i], &row[cj ^ i], &mone);
      }

#else

      if (scal == 0) {
        field_add(&row[cj ^ i], &row[cj ^ i], &one);
      } else if (scal == 1) {
        field_add(&row[cj ^ i], &row[cj ^ i], &onexi);
      } else if (scal == 2) {
        field_add(&row[cj ^ i], &row[cj ^ i], &mone);
      } else {
        field_add(&row[cj ^ i], &row[cj ^ i], &monexi);
      }

#endif
    }

    for (int j = 0; j < 16; j++) {
      if (!field_is_zero(&row[j])) {
        return i;
      }
    }
  }
}

void compute_change_theta_coords_dim4(
    field_mat_16x16_t *mat, const mod4_mat_4x4_t *A, const mod4_mat_4x4_t *B,
    const mod4_mat_4x4_t *C, const mod4_mat_4x4_t *D, const unsigned int is_neg,
    const int nv_ind) {
  unsigned int i0;
  if (nv_ind == -1) {
    // Case of unkown non vanishing index. Not constant time.
    field_t row_i0[16];
    i0 = choose_non_vanishing_index(row_i0, C, D, is_neg);
    for (int j = 0; j < 16; j++) {
      field_copy(&(*mat)[i0][j], &row_i0[j]);
    }
  } else {
    // Case of kown non vanishing index. Constant time guarantee.
    i0 = nv_ind;
  }

#ifdef FP_ONLY
  field_t one, mone;
  field_set_one(&one);
  field_neg(&mone, &one);
  unsigned int ctrl;
#else
  field_t one, onexi, monexi;
  field_set_one(&one);
  fp_set_zero(&onexi.re);
  fp_set_one(&onexi.im);
  fp_set_zero(&monexi.re);
  fp_neg(&monexi.im, &onexi.im);
  field_cswap(&onexi, &monexi, is_neg);
  unsigned int ctrl1, ctrl2;
  field_t neg_sel;
#endif

  mod2_vec_4_t I0, I, J, AI, BI, AIpCJ, BIpDJ;
  index_to_multindex_dim4(I0, i0);
  unsigned int scal, aipcjpi0;
  field_t sel;
  for (int i = 0; i < 16; i++) {
    if ((nv_ind != -1) || (i != i0)) {
      index_to_multindex_dim4(I, i);
      mod4_mat_4x4_eval(&AI, A, &I);
      mod4_mat_4x4_eval(&BI, B, &I);

      for (int j = 0; j < 16; j++) {
        field_set_zero(&(*mat)[i][j]);
      }

      for (int j = 0; j < 16; j++) {
        index_to_multindex_dim4(J, j);
        mod4_mat_4x4_eval(&AIpCJ, C, &J);
        mod4_mat_4x4_eval(&BIpDJ, D, &J);
        mod4_add_4(&AIpCJ, &AIpCJ, &AI);
        mod4_add_4(&BIpDJ, &BIpDJ, &BI);
        scal = (mod4_scal_prod_4(I, J) - mod4_scal_prod_4(&AIpCJ, &BIpDJ) -
                2 * mod4_scal_prod_4(&I0, &BIpDJ)) & 3;

        aipcjpi0 = multindex_to_index_dim4(AIpCJ) ^ i0;

#ifdef FP_ONLY
        ctrl = (scal == 2);
        field_select(&sel, &one, &mone, ctrl);
#else
        ctrl1 = (scal & 1 == 1);  // Choose i = onexi
        ctrl2 = (!ctrl1 & (scal == 2)) || (ctrl1 & (scal == 3));  // Select the
        // negative option
        field_select(&sel, &one, &onexi, ctrl1);
        field_neg(&neg_sel, &sel);
        field_select(&sel, &sel, &neg_sel, ctrl2);
#endif

        field_add(&(*mat)[i][aipcjpi0], &(*mat)[i][aipcjpi0], &sel);
      }
    }
  }
}

void couple_theta_point_dim2_to_theta_point_dim4(
    theta_point_dim4_t *res, const couple_theta_point_dim2_t *in,
    const field_mat_16x16_t *mat_change_theta_coords) {
  field_t product[16], tmp;
  field_mul(&product[0], &in->P1.x, &in->P2.x);
  field_mul(&product[1], &in->P1.y, &in->P2.x);
  field_mul(&product[2], &in->P1.z, &in->P2.x);
  field_mul(&product[3], &in->P1.t, &in->P2.x);
  field_mul(&product[4], &in->P1.x, &in->P2.y);
  field_mul(&product[5], &in->P1.y, &in->P2.y);
  field_mul(&product[6], &in->P1.z, &in->P2.y);
  field_mul(&product[7], &in->P1.t, &in->P2.y);
  field_mul(&product[8], &in->P1.x, &in->P2.z);
  field_mul(&product[9], &in->P1.y, &in->P2.z);
  field_mul(&product[10], &in->P1.z, &in->P2.z);
  field_mul(&product[11], &in->P1.t, &in->P2.z);
  field_mul(&product[12], &in->P1.x, &in->P2.t);
  field_mul(&product[13], &in->P1.y, &in->P2.t);
  field_mul(&product[14], &in->P1.z, &in->P2.t);
  field_mul(&product[15], &in->P1.t, &in->P2.t);

  for (int i = 0; i < 16; i++) {
    field_set_zero(&res->coords[i]);

    for (int j = 0; j < 16; j++) {
      field_mul(&tmp, &product[j], &(*mat_change_theta_coords)[i][j]);
      field_add(&res->coords[i], &tmp);
    }
  }
}