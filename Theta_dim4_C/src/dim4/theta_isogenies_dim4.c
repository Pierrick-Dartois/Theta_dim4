#include <theta_isogenies_dim4.h>
#include <theta_structures_dim4.h>

void isogeny_eval_dim4(theta_point_dim4_t *image, const theta_point_dim4_t *point, const theta_struct_dim4_t *codomain){
	theta_dim4_sq(image->coords,point->coords);
	hadamard_dim4(image->coords,image->coords);
	//for(int i=0; i<16; i++){
		//field_mul(image->coords);
	//}
}

