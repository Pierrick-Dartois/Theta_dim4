#include <theta_structures_dim4.h>

static inline void theta_copy(theta_point_dim4_t *out, const theta_point_dim4_t *in){
	field_copy(&out->coords[i],&in->coords[i]);
}

static inline void hadamard_2(field_t *out, const field_t *in) {
  field_add(&out[0], &in[0], &in[1]);
  field_sub(&out[1], &in[0], &in[1]);
}

static inline void hadamard_4(field_t *out, const field_t *in) {
  field_t in1[2], in2[2], out1[2], out2[2];

  field_copy(&in1[0], &in[0]);
  field_copy(&in1[1], &in[1]);
  field_copy(&in2[0], &in[2]);
  field_copy(&in2[1], &in[3]);

  hadamard_2(out1, in1);
  hadamard_2(out2, in2);

  field_add(&out[0], &out1[0], &out2[0]);
  field_add(&out[1], &out1[1], &out2[1]);
  field_sub(&out[2], &out1[0], &out2[0]);
  field_sub(&out[3], &out1[1], &out2[1]);
}

static inline void hadamard_8(field_t *out, const field_t *in) {
  field_t in1[4], in2[4], out1[4], out2[4];

  for (int i = 0; i < 4; i++) {
    field_copy(&in1[i], &in[i]);
    field_copy(&in2[i], &in[i + 4]);
  }

  hadamard_4(out1, in1);
  hadamard_4(out2, in2);

  for (int i = 0; i < 4; i++) {
    field_add(&out[i], &out1[i], &out2[i]);
    field_sub(&out[i + 4], &out1[i], &out2[i]);
  }
}

static inline void hadamard_16(field_t *out, const field_t *in) {
  field_t in1[8], in2[8], out1[8], out2[8];

  for (int i = 0; i < 8; i++) {
    field_copy(&in1[i], &in[i]);
    field_copy(&in2[i], &in[i + 8]);
  }

  hadamard_4(out1, in1);
  hadamard_4(out2, in2);

  for (int i = 0; i < 8; i++) {
    field_add(&out[i], &out1[i], &out2[i]);
    field_sub(&out[i + 8], &out1[i], &out2[i]);
  }
}

void hadamard_dim4(field_t *out, const field_t *in) { hadamard_16(out, in); }

void theta_dim4_sq(field_t *out, const field_t *in) {
  for (int i = 0; i < 16; i++) {
    field_sqr(&out[i], &in[i]);
  }
}

void theta_struct_arith_precomp(theta_struct_dim4_t *theta_struct) {
  for (int i = 0; i < 16; i++) {
    field_copy(&theta_struct->null_point.coords[i],
    	&theta_struct->inv_null_point[i]);
  }
  field_proj_batched_inv(theta_struct->null_point.coords, 16);

  theta_dim4_sq(theta_struct->inv_codomain_dual_null_point_sq,
                theta_struct->null_point.coords);
  hadamard_dim4(theta_struct->inv_codomain_dual_null_point_sq,
                theta_struct->inv_codomain_dual_null_point_sq);
  field_proj_batched_inv(theta_struct->inv_codomain_dual_null_point_sq, 16);

  theta_struct->arith_precomp = 1;
}

void theta_double(theta_point_dim4_t *out, const theta_point_dim4_t *in,
                  theta_struct_dim4_t *theta_struct) {
	if (!theta_struct->arith_precomp){
		theta_struct_arith_precomp(theta_struct);
	}

	theta_dim4_sq(out->coords,in->coords);
	hadamard_dim4(out->coords,out->coords);
	theta_dim4_sq(out->coords,out->coords);

	for(int i=0; i<16; i++) {
		field_mul(&out->coords[i],&out->coords[i],&theta_struct->inv_codomain_dual_null_point_sq[i]);
	}
	
	hadamard_dim4(out->coords,out->coords);

	for(int i=0; i<16; i++) {
		field_mul(&out->coords[i],&out->coords[i],&theta_struct->inv_null_point[i]);
	}
}

void theta_double_iter(theta_point *out, const theta_point_dim4_t *in,
                  theta_struct_dim4_t *theta_struct, int n) {

	theta_copy(out, in);
	for(int i=0; i<n; i++){
		theta_double(out, out,theta_struct);
	}
}



