#include "isog.h"
#include <field.h>
#include <assert.h>

// since we use degree 4 isogeny steps, we need to handle the odd case with care
static uint32_t
ec_eval_even_strategy(ec_curve_t *curve,
                      ec_point_t *points,
                      unsigned len_points,
                      const ec_point_t *kernel,
                      const int isog_len)
{
    ec_curve_normalize_A24(curve);
    ec_point_t A24;
    copy_point(&A24, &curve->A24);

    int space = 1;
    for (int i = 1; i < isog_len; i *= 2)
        ++space;

    // Stack of remaining kernel points and their associated orders
    ec_point_t splits[space];
    uint16_t todo[space];
    splits[0] = *kernel;
    todo[0] = isog_len;

    int current = 0; // Pointer to current top of stack

    // Chain of 4-isogenies
    for (int j = 0; j < isog_len / 2; ++j) {
        assert(current >= 0);
        assert(todo[current] >= 1);
        // Get the next point of order 4
        while (todo[current] != 2) {
            assert(todo[current] >= 3);
            // A new split will be added
            ++current;
            assert(current < space);
            // We set the seed of the new split to be computed and saved
            copy_point(&splits[current], &splits[current - 1]);
            // if we copied from the very first element, then we perform one additional doubling
            unsigned num_dbls = todo[current - 1] / 4 * 2 + todo[current - 1] % 2;
            todo[current] = todo[current - 1] - num_dbls;
            while (num_dbls--)
                xDBL_A24(&splits[current], &splits[current], &A24, false);
        }

        if (j == 0) {
            assert(field_is_one(&A24.z));
            if (!ec_is_four_torsion(&splits[current], curve))
                return -1;

            ec_point_t T;
            xDBL_A24(&T, &splits[current], &A24, false);
            if (field_is_zero(&T.x))
                return -1; // special isogenies not allowed
        } else {
            assert(todo[current] == 2);
#ifndef NDEBUG
            if (field_is_zero(&splits[current].z))
                debug_print("splitting point z coordinate is unexpectedly zero");

            ec_point_t test;
            xDBL_A24(&test, &splits[current], &A24, false);
            if (field_is_zero(&test.z))
                debug_print("z coordinate is unexpectedly zero before doubling");
            xDBL_A24(&test, &test, &A24, false);
            if (!field_is_zero(&test.z))
                debug_print("z coordinate is unexpectedly not zero after doubling");
#endif
        }

        // Evaluate 4-isogeny
        ec_kps4_t kps4;
        xisog_4(&kps4, &A24, splits[current]);
        xeval_4(splits, splits, current, &kps4);
        for (int i = 0; i < current; ++i)
            todo[i] -= 2;
        xeval_4(points, points, len_points, &kps4);

        --current;
    }
    assert(isog_len % 2 ? !current : current == -1);

    // Final 2-isogeny
    if (isog_len % 2) {
#ifndef NDEBUG
        if (field_is_zero(&splits[0].z))
            debug_print("splitting point z coordinate is unexpectedly zero");
        ec_point_t test;
        copy_point(&test, &splits[0]);
        xDBL_A24(&test, &test, &A24, false);
        if (!field_is_zero(&test.z))
            debug_print("z coordinate is unexpectedly not zero after doubling");
#endif

        // We need to check the order of this point in case there were no 4-isogenies
        if (isog_len == 1 && !ec_is_two_torsion(&splits[0], curve))
            return -1;
        if (field_is_zero(&splits[0].x)) {
            // special isogenies not allowed
            // this case can only happen if isog_len == 1; otherwise the
            // previous 4-isogenies we computed ensure that $T=(0:1)$ is put
            // as the kernel of the dual isogeny
            return -1;
        }

        ec_kps2_t kps2;
        xisog_2(&kps2, &A24, splits[0]);
        xeval_2(points, points, len_points, &kps2);
    }

    // Output curve in the form (A:C)
    A24_to_AC(curve, &A24);

    curve->is_A24_computed_and_normalized = false;

    return 0;
}

uint32_t
ec_eval_even(ec_curve_t *image, ec_isog_even_t *phi, ec_point_t *points, unsigned len_points)
{
    copy_curve(image, &phi->curve);
    return ec_eval_even_strategy(image, points, len_points, &phi->kernel, phi->length);
}

// naive implementation
uint32_t
ec_eval_small_chain(ec_curve_t *curve,
                    const ec_point_t *kernel,
                    int len,
                    ec_point_t *points,
                    unsigned len_points,
                    bool special) // do we allow special isogenies?
{

    ec_point_t A24;
    AC_to_A24(&A24, curve);

    ec_kps2_t kps;
    ec_point_t small_K, big_K;
    copy_point(&big_K, kernel);

    for (int i = 0; i < len; i++) {
        copy_point(&small_K, &big_K);
        // small_K = big_K;
        for (int j = 0; j < len - i - 1; j++) {
            xDBL_A24(&small_K, &small_K, &A24, false);
        }
        // Check the order of the point before the first isogeny step
        if (i == 0 && !ec_is_two_torsion(&small_K, curve))
            return (uint32_t)-1;
        // Perform isogeny step
        if (field_is_zero(&small_K.x)) {
            if (special) {
                ec_point_t B24;
                xisog_2_singular(&kps, &B24, A24);
                xeval_2_singular(&big_K, &big_K, 1, &kps);
                xeval_2_singular(points, points, len_points, &kps);
                copy_point(&A24, &B24);
            } else {
                return (uint32_t)-1;
            }
        } else {
            xisog_2(&kps, &A24, small_K);
            xeval_2(&big_K, &big_K, 1, &kps);
            xeval_2(points, points, len_points, &kps);
        }
    }
    A24_to_AC(curve, &A24);

    curve->is_A24_computed_and_normalized = false;
    return 0;
}

uint32_t
ec_isomorphism(ec_isom_t *isom, const ec_curve_t *from, const ec_curve_t *to)
{
    field_t t0, t1, t2, t3, t4;

    field_mul(&t0, &from->A, &from->C);
    field_mul(&t1, &to->A, &to->C);

    field_mul(&t2, &t1, &to->C); // toA*toC^2
    field_add(&t3, &t2, &t2);
    field_add(&t3, &t3, &t3);
    field_add(&t3, &t3, &t3);
    field_add(&t2, &t2, &t3); // 9*toA*toC^2
    field_sqr(&t3, &to->A);
    field_mul(&t3, &t3, &to->A); // toA^3
    field_add(&t3, &t3, &t3);
    field_sub(&isom->Nx, &t3, &t2); // 2*toA^3-9*toA*toC^2
    field_mul(&t2, &t0, &from->A);  // fromA^2*fromC
    field_sqr(&t3, &from->C);
    field_mul(&t3, &t3, &from->C); // fromC^3
    field_add(&t4, &t3, &t3);
    field_add(&t3, &t4, &t3);             // 3*fromC^3
    field_sub(&t3, &t3, &t2);             // 3*fromC^3-fromA^2*fromC
    field_mul(&isom->Nx, &isom->Nx, &t3); // lambda_x = (2*toA^3-9*toA*toC^2)*(3*fromC^3-fromA^2*fromC)

    field_mul(&t2, &t0, &from->C); // fromA*fromC^2
    field_add(&t3, &t2, &t2);
    field_add(&t3, &t3, &t3);
    field_add(&t3, &t3, &t3);
    field_add(&t2, &t2, &t3); // 9*fromA*fromC^2
    field_sqr(&t3, &from->A);
    field_mul(&t3, &t3, &from->A); // fromA^3
    field_add(&t3, &t3, &t3);
    field_sub(&isom->D, &t3, &t2); // 2*fromA^3-9*fromA*fromC^2
    field_mul(&t2, &t1, &to->A);   // toA^2*toC
    field_sqr(&t3, &to->C);
    field_mul(&t3, &t3, &to->C); // toC^3
    field_add(&t4, &t3, &t3);
    field_add(&t3, &t4, &t3);           // 3*toC^3
    field_sub(&t3, &t3, &t2);           // 3*toC^3-toA^2*toC
    field_mul(&isom->D, &isom->D, &t3); // lambda_z = (2*fromA^3-9*fromA*fromC^2)*(3*toC^3-toA^2*toC)

    // Mont -> SW -> SW -> Mont
    field_mul(&t0, &to->C, &from->A);
    field_mul(&t0, &t0, &isom->Nx); // lambda_x*toC*fromA
    field_mul(&t1, &from->C, &to->A);
    field_mul(&t1, &t1, &isom->D);  // lambda_z*fromC*toA
    field_sub(&isom->Nz, &t0, &t1); // lambda_x*toC*fromA - lambda_z*fromC*toA
    field_mul(&t0, &from->C, &to->C);
    field_add(&t1, &t0, &t0);
    field_add(&t0, &t0, &t1);             // 3*fromC*toC
    field_mul(&isom->D, &isom->D, &t0);   // 3*lambda_z*fromC*toC
    field_mul(&isom->Nx, &isom->Nx, &t0); // 3*lambda_x*fromC*toC

    return (field_is_zero(&isom->Nx) | field_is_zero(&isom->D));
}

void
ec_iso_eval(ec_point_t *P, ec_isom_t *isom)
{
    field_t tmp;
    field_mul(&P->x, &P->x, &isom->Nx);
    field_mul(&tmp, &P->z, &isom->Nz);
    field_add(&P->x, &P->x, &tmp);
    field_mul(&P->z, &P->z, &isom->D);
}

void
ec_2_isog_chain(ec_2_isog_chain_t *chain,const ec_point_t *kernel, const ec_curve_t *domain, unsigned int len,
    const unsigned int *strategy)
{
    ec_kps2_t kps[len];
    ec_point_t A24[len+1];
    bool is_singular[len];

    unsigned int tmp,
        log2_of_e, // Height of the strategy tree topology
        strat_idx = 0, // Current element of the strategy to be used
        block = 0, // Keeps track of point order
        current = 0; // Number of points being carried

    for (tmp = len, log2_of_e = 0; tmp > 0; tmp >>= 1, ++log2_of_e)
        ;
    log2_of_e *= 2;

    ec_point_t kernel_elements[log2_of_e];
    unsigned int XDBLs[log2_of_e]; // Number of doubles performed

    AC_to_A24(&A24[0],domain);
    copy_point(&kernel_elements[0], kernel);

    for(int k=0;k<len;k++){
        while(block!=len-1-k){
            current += 1;

            // Append the last kernel element and performs the doublings
            copy_point(&kernel_elements[current],&kernel_elements[current-1]);
            for(int j=0;j<strategy[strat_idx];j++){
                xDBL_A24(&kernel_elements[current],&kernel_elements[current],&A24[k]);
            }

            // Update bookkeeping variables
            XDBLs[current]=strategy[strat_idx];
            block+=strategy[strat_idx];
            strat_idx+=1;
        }

        if(field_is_zero(&kernel_elements[current].x)){
            xisog_2_singular(&kps[k], &A24[k+1], A24[k]);
            xeval_2_singular(kernel_elements, kernel_elements, current, &kps[k]);
            is_singular[k]=true;
        }
        else{
            xisog_2(&kps[k], &A24[k+1], kernel_elements[current]);
            xeval_2(kernel_elements, kernel_elements, current, &kps[k]);
            is_singular[k]=false;
        }

        block -= XDBLs[current];
        XDBLs[current] = 0;
        current -= 1;
    }

    chain->len=len;
    chain->kps=(ec_kps2_t *)malloc(len*sizeof(ec_kps2_t));
    chain->is_singular=(bool *)malloc(len*sizeof(bool));
    chain->A24=(ec_point_t *)malloc((len+1)*sizeof(ec_point_t));
    //kps;
    for(int i=0;i<len;i++){
        copy_point(&chain->kps[i].K,&kps[i].K);
        chain->is_singular[i]=is_singular[i];
        copy_point(&chain->A24[i],&A24[i]);
    }
    copy_point(&chain->A24[len],&A24[len]);
    copy_curve(&chain->domain,domain);
    ec_curve_init(&chain->codomain);
    A24_to_AC(&chain->codomain,&A24[len]);
}

void 
ec_eval_2_isog_chain(ec_point_t *Q, const ec_point_t *P, const ec_2_isog_chain_t *chain)
{
    copy_point(Q,P);
    for(int i=0;i<chain->len;i++){
        if(chain->is_singular[i]){
            xeval_2_singular(Q, Q, 1, &chain->kps[i]);
        }
        else{
            xeval_2(Q, Q, 1, &chain->kps[i]);
        }
    }
}

void
del_2_isog_chain(ec_2_isog_chain_t *chain){
    free(chain->kps);
    free(chain->is_singular);
    free(chain->A24);
}

void
ec_2_torsion_point(ec_point_t *P, const ec_curve_t *E)
{
    // Returns a non-zero 2-torsion point of E !=(0,0).
    field_t t0, t1, t2;

    field_add(&t0,&E->C,&E->C); // 2C
    field_add(&t1,&E->A,&t0); // A+2C
    field_sub(&t2,&E->A,&t0); // A-2C
    field_mul(&t1,&t1,&t2); // (A+2C)(A-2C)
    field_sqrt(&t1); // \sqrt((A+2C)(A-2C))
    field_add(&P->x,&t1,&E->A); // x = A+\sqrt((A+2C)(A-2C))
    field_neg(&P->z,&t0); // z = -2C
}

void
ec_odd_isog_chain(ec_odd_isog_chain_t *chain,const ec_point_t *kernel, const ec_curve_t *domain, unsigned int l, 
    unsigned int len, const unsigned int *strategy)
{
    ec_kps_t kps[len];
    ec_point_t A24[len+1], A3, P2;
    digit_t tabl[1];
    int nbits_l=nbits_int(l), d=(l-1)/2;

    if(l!=3){
        ec_2_torsion_point(&P2,domain);
        // A point of 2-torsion should be propagated through the chain when l>3
        mp_set_small(tabl,l,1);
    }

    

    unsigned int tmp,
        log2_of_e, // Height of the strategy tree topology
        strat_idx = 0, // Current element of the strategy to be used
        block = 0, // Keeps track of point order
        current = 0; // Number of points being carried

    for (tmp = len, log2_of_e = 0; tmp > 0; tmp >>= 1, ++log2_of_e)
        ;
    log2_of_e *= 2;

    ec_point_t kernel_elements[log2_of_e];
    unsigned int XDBLs[log2_of_e]; // Number of multiplications performed

    AC_to_A24(&A24[0],domain);
    copy_point(&kernel_elements[0], kernel);

    for(int k=0;k<len;k++){
        while(block!=len-1-k){
            current += 1;

            // Append the last kernel element and performs the l-multiplications
            copy_point(&kernel_elements[current],&kernel_elements[current-1]);
            if(l==3){
                field_copy(&A3.x,&A24[k].x);// A+2C
                field_sub(&A3.z,&A24[k].x,&A24[k].z);// A-2C = A+2C - 4C
                for(int j=0;j<strategy[strat_idx];j++){
                    xTPL(&kernel_elements[current],&kernel_elements[current],&A3);
                }
            }
            else{
                for(int j=0;j<strategy[strat_idx];j++){
                    xMUL_A24(&kernel_elements[current],&kernel_elements[current],tabl,nbits_l,&A24[k]);
                }
            }

            // Update bookkeeping variables
            XDBLs[current]=strategy[strat_idx];
            block+=strategy[strat_idx];
            strat_idx+=1;
        }



        if(l==3){
            xisog_3(&kps[k], &A24[k+1], kernel_elements[current]);
            xeval_3(kernel_elements, kernel_elements, current, &kps[k]);
        }
        else{
            xisog_odd(&kps[k], &A24[k+1], kernel_elements[current], &A24[k], P2, d);
            xeval_odd(&P2, &P2, 1, &kps[k], d);
            xeval_odd(kernel_elements, kernel_elements, current, &kps[k], d);
        }

        block -= XDBLs[current];
        XDBLs[current] = 0;
        current -= 1;
    }

    chain->len=len;
    chain->d=d;
    chain->kps=(ec_kps_t *)malloc(len*sizeof(ec_kps_t));
    chain->A24=(ec_point_t *)malloc((len+1)*sizeof(ec_point_t));
    //kps;
    for(int i=0;i<len;i++){
        chain->kps[i]=kps[i];
        copy_point(&chain->A24[i],&A24[i]);
    }
    copy_point(&chain->A24[len],&A24[len]);
    copy_curve(&chain->domain,domain);
    ec_curve_init(&chain->codomain);
    A24_to_AC(&chain->codomain,&A24[len]);
}

void 
ec_eval_odd_isog_chain(ec_point_t *Q, const ec_point_t *P, const ec_odd_isog_chain_t *chain)
{
    copy_point(Q,P);
    for(int i=0;i<chain->len;i++){
        if(chain->d==1){
            xeval_3(Q, Q, 1, &chain->kps[i]);
        }
        else{
            xeval_odd(Q, Q, 1, &chain->kps[i], chain->d);
        }
    }
}

void
del_odd_isog_chain(ec_odd_isog_chain_t *chain){
    for(int i=0;i<chain->len;i++){
        free(chain->kps[i].K);
    }
    free(chain->kps);
    free(chain->A24);
}
