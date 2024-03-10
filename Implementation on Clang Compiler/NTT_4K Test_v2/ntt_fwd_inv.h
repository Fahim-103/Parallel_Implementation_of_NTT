#include "ntt_calc.h"

//
// Created by secure on 11/7/23.
//

#ifndef WORD_SIZE
#define WORD_SIZE 16 // "defs.h" x >> 64
#endif

#ifndef NTT_FWD_INV_H
#define NTT_FWD_INV_H
#endif //NTT_FWD_INV_H


#define HIGH_WORD(x) ((x) >> WORD_SIZE)
#define LOW_WORD(x)  ((x)&WORD_SIZE_MASK)

#ifndef L_HIGH_WORD
#  define L_HIGH_WORD HIGH_WORD // "defs.h" x >> 64
#endif

////////// function required for havybtrfly///////////////////

////// 4q-2q-q //////////// modulus optimization//////////////

static inline uint64_t reduce_4q_to_2q(const uint64_t val, const uint64_t q)
{
    return (val < 2 * q) ? val : val - 2 * q;  // value < q is a must, at random_buf(), values[i] = rand() % q
}

/////////////////////////////////////////////////////////////////////
static inline uint64_t reduce_2q_to_q(const uint64_t val, const uint64_t q)
{
    return (val < q) ? val : val - q;
}

static inline uint64_t reduce_4q_to_q(const uint64_t val, const uint64_t q)
{
    return reduce_2q_to_q(reduce_4q_to_2q(val, q), q);
}
////////////////////////////////////////////////////////////////////////////////

///////////////////fast_mul_mod_q2/////////////////////////////////////////

static inline uint64_t fast_mul_mod_q2(const mul_op_t w, const uint64_t t, const uint64_t q) // w_power[m+i] is in the mul_op_t str, a_original[i + (N/2)], q is the prime //generates shoupmodmul()
{                                                                     // only X1-X2 is not passed, direct X2 is passed. Negation is is done in harvey_fwd_butterfly() last line
    const uint64_t Q = L_HIGH_WORD(w.con * t); //w.con 128 bit, t=a_original[i + (N/2)] is 64 bit, first 64 bits are taken// w.con = w' = (wB/q)
    return w.op * t - Q * q; // t = Y, (w.op*Y-w.con*Y*q) = w*Y mod q, has been performed
}

///////////////////fast_mul_mod_q - for calculating inverse  /////////////////////////////////////////


/////////////////////////////////////////////////////////////


static inline void
harvey_fwd_butterfly(uint64_t *X, uint64_t *Y, const mul_op_t w, const uint64_t q)
{
    const uint64_t q2 = q << 1;
    const uint64_t X1 = reduce_4q_to_2q(*X + *Y, q);
    const uint64_t T  = *X - *Y + q2;

    *X = X1;
    *Y = fast_mul_mod_q2(w, T, q);
}



void fwd_ntt_ref_harvey_lazy(uint64_t       a[], // randomly generated x(n), number of elements in a[] = N
                             const uint64_t N,   // N = 2^m, size of x(n)
                             const uint64_t q,   // Zq or Rq, ring of q. q is a prime
                             const uint64_t w[], // t->w_powers.ptr, omega,w_powers ?(NbyN) 2D matrix?, converted to 1D array, pwrs of w // created at "test_cases.h" as a part test_case_t
                             const uint64_t w_con[])  //w converted
{
    size_t t = N;

    for(size_t m = 1; m < N; m <<= 1) // m = m*2, m=m<<1, t = t/2; number of blocks
    {
        t = t >> 1;
        for(size_t i = 0; i < m; i++) {       //block-wise iteration

               //w.con is set here // (t->w_powers_con.ptr) from fwd_ntt_ref_harvey() from test_correctness() from "test_correctness.c"
            int k= 0;

           // LOOP_UNROLL_4  //4 times unroll the loop for faster execution
            for(size_t j = 2*i*t; j < (2*i*t + t); j++) {
               // printf("block created with a[%zu] and a[%ld] and w_power[%ld]\n", j, j+t, m+i);

                const mul_op_t w1 = {w[m*k], w_con[m*k]};
               // printf("from fwd_ntt_ref_harvey_lazy - taking a[%zu] and a[%zu] with w1[%zu] of value %lu\n",j, j+t, k, w1.op);

              //  printf("--> from fwd_ntt_ref_harvey_lazy - taking a[%zu]=%ld and a[%zu]=%ld with w1[%zu] of value %lu\n",j,a[j], j+t ,a[j+t], m*k, w1.op);
                harvey_fwd_butterfly(&a[j], &a[j + t], w1, q);  //radix-2 value of a[j], a[j+N/2] calculated,  function defn at "fast_mul_operators.h"
                //printf("<-- from fwd_ntt_ref_harvey_lazy - generating a[%zu]=%ld and a[%zu]=%ld with w1[%zu] of value %lu\n",j,a[j], j+t ,a[j+t], m*k, w1.op);

                k = k+1;
            }
          //  printf("loop ends for i = %zu and m= %ld\n\n\n", i, m);
        }
    }
}


static inline void fwd_ntt_ref_harvey(uint64_t       a[],
                                      const uint64_t N,
                                      const uint64_t q,
                                      const uint64_t w[],
                                      const uint64_t w_con[],
                                      const uint64_t width)
{
    fwd_ntt_ref_harvey_lazy(a, N, q, w, w_con);

    // Final reduction
    uint64_t temp[N];
    for(size_t i = 0; i < N; i++) {

       // printf ("input a[%zu] %lu ", i, a[i]);

        temp[i] = reduce_4q_to_q(a[i], q);

       // printf ("output a[%zu] %lu ****** ", i, a[i]);
      //  printf("inside lazy function: %ld\n", a[i]);

    }

    for(size_t i = 0; i < N; i++) {
        a[bit_rev_idx(i, width)] = temp[i];
    }

   // printf("\n");

}



////////////////////////////Inverse NTT Functions///////////////////////////////////////////////////////////

////////////////////////supporting funcitons///////////////////////////////////////////////////

static inline uint64_t
fast_mul_mod_q(const mul_op_t w, const uint64_t t, const uint64_t q) // w_powers, a_original[i + (N/2)], q is the prime
{
    return reduce_2q_to_q(fast_mul_mod_q2(w, t, q), q); //required to find (1/N)*[]
}

/////////////////// multiplying with 1/N ////////////////////////////////////

static inline void harvey_bkw_butterfly_final(uint64_t *     X, // a[j] comes from bkw_butterfly
                                              const mul_op_t n_inv,
                                              const uint64_t q)
{

    *X = fast_mul_mod_q2(n_inv, *X, q);   //instead of dividing by N, we do modulus operation

}


/////////function called in main//////////////////////////

void inv_ntt_ref_harvey(uint64_t       a[], // after fwd_ntt a[]
                        const uint64_t N, // N = 2^m, size of x(n)
                        const uint64_t q, // prime
                        const mul_op_t n_inv, // manual input n_inv.op found, n_inv.con not found
                        const uint64_t word_size, // 64
                        const uint64_t w[], // t->w_inv_powers.ptr
                        const uint64_t w_con[],
                        const uint64_t width)

{

    size_t t = N;

    for(size_t m = 1; m < N; m <<= 1) // m = m*2, m=m<<1, t = t/2; number of blocks
    {
        t = t >> 1;
        for(size_t i = 0; i < m; i++) {       //block-wise iteration

            //w.con is set here // (t->w_powers_con.ptr) from fwd_ntt_ref_harvey() from test_correctness() from "test_correctness.c"
            int k= 0;

            // LOOP_UNROLL_4  //4 times unroll the loop for faster execution
            for(size_t j = 2*i*t; j < (2*i*t + t); j++) {
                // printf("block created with a[%zu] and a[%ld] and w_power[%ld]\n", j, j+t, m+i);

                const mul_op_t w1 = {w[m*k], w_con[m*k]};
               // printf("from bkwd_ntt_ref_harvey_lazy - taking a[%zu]=%ld and a[%zu]=%ld with w1[%zu] of value %lu\n",j,a[j], j+t ,a[j+t], m*k, w1.op);
                harvey_fwd_butterfly(&a[j], &a[j + t], w1, q);  //radix-2 value of a[j], a[j+N/2] calculated,  function defn at "fast_mul_operators.h"
               // printf("from bkwd_ntt_ref_harvey_lazy - generating a[%zu]=%ld and a[%zu]=%ld with w1[%zu] of value %lu\n",j,a[j], j+t ,a[j+t], m*k, w1.op);
                k = k+1;
            }

        }
    }

    uint64_t temp_inv[N];

    for(size_t j = 0; j < N; j++) {
        harvey_bkw_butterfly_final(&a[j], n_inv, q); // revised value of a[j], a[j+N/2] calculated,   function defn at "fast_mul_operators.h"
        //temp_inv[j] = reduce_4q_to_q(a[j], q);
        temp_inv[j] = a[j] % q;
    }

    for(size_t k=0; k< N; k++){

        a[bit_rev_idx(k, width)] = temp_inv[k];

    }

}






