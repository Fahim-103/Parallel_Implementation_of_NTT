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

//fwd_ntt_ref_harvey(a, t->n, t->q, t->w_powers.ptr, t->w_powers_con.ptr);////fwd ntt call

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
    return w.op * t - Q * q; // T at harvey_fwd_butterfly not fwd_ref.........
}

///////////////////fast_mul_mod_q - for calculating inverse  /////////////////////////////////////////


/////////////////////////////////////////////////////////////


static inline void
harvey_fwd_butterfly(uint64_t *X, uint64_t *Y, const mul_op_t w, const uint64_t q) // a_original[i] , a_original[i+N/2]
// w is a struct with w.op, w.con
{
    const uint64_t q2 = q << 1;  //q2 = 2q
    const uint64_t X1 = reduce_4q_to_2q(*X, q); // input a[i]%q... a[1] = 107 for first random generation, q 251
    // X1 = *X as per random_buff generation of value
    const uint64_t T  = fast_mul_mod_q2(w, *Y, q); // w = w_powers[m+i], *Y = a_orig[i+N/2]

    *X = X1 + T;      // implementation of Harvey Butt, a[i] value changes to fwd_ntt
    *Y = X1 - T + q2; // implementation of Harvey Butt, a[i+N/2]
}



void fwd_ntt_ref_harvey_lazy(uint64_t       a[], // randomly generated x(n), number of elements in a[] = N
                             const uint64_t N,   // N = 2^m, size of x(n)
                             const uint64_t q,   // Zq or Rq, ring of q. q is a prime
                             const uint64_t w[], // t->w_powers.ptr, omega,w_powers ?(NbyN) 2D matrix?, converted to 1D array, pwrs of w // created at "test_cases.h" as a part test_case_t
                             const uint64_t w_con[])  //w'
{
    size_t t = N >> 1; // t = N/2

    for(size_t m = 1; m < N; m <<= 1, t >>= 1) // m = m*2, m=m<<1, t = t/2
    {
        size_t k = 0;
        for(size_t i = 0; i < m; i++) {
            const mul_op_t w1 = {w[m + i], w_con[m + i]};   //w.con is set here // (t->w_powers_con.ptr) from fwd_ntt_ref_harvey() from test_correctness() from "test_correctness.c"

            // LOOP_UNROLL_4  //4 times unroll the loop for faster execution
            for(size_t j = k; j < k + t; j++) {
                harvey_fwd_butterfly(&a[j], &a[j + t], w1, q);  //radix-2 value of a[j], a[j+N/2] calculated,  function defn at "fast_mul_operators.h"
            }
            k = k + (2 * t);
        }
    }
} // lazy_version + reduction = actual_function // called at "test_corectness.c in test_radix2_scalar()"


static inline void fwd_ntt_ref_harvey(uint64_t       a[],
                                      const uint64_t N,
                                      const uint64_t q,
                                      const uint64_t w[],
                                      const uint64_t w_con[])
{
    fwd_ntt_ref_harvey_lazy(a, N, q, w, w_con);

    // Final reduction
    for(size_t i = 0; i < N; i++) {
        a[i] = reduce_4q_to_q(a[i], q);  //from "fast_mul_operators.h" all a[i] are q modulused

        //  printf("inside lazy function: %ld\n", a[i]);

    }
}

////////////////////////////Inverse NTT Functions///////////////////////////////////////////////////////////

////////////////////////supporting funcitons///////////////////////////////////////////////////

static inline uint64_t
fast_mul_mod_q(const mul_op_t w, const uint64_t t, const uint64_t q) // w_powers, a_original[i + (N/2)], q is the prime
{
    return reduce_2q_to_q(fast_mul_mod_q2(w, t, q), q); //required to find (1/N)*[]
}

/////////////////////////////////////////////////////////////////////////////////////////////

/////////radix-2 bkwd butterfly//////////

static inline void
harvey_bkw_butterfly(uint64_t *X, uint64_t *Y, const mul_op_t w, const uint64_t q)
{
    const uint64_t q2 = q << 1;
    const uint64_t X1 = reduce_4q_to_2q(*X + *Y, q);
    const uint64_t T  = *X - *Y + q2;

    *X = X1;
    *Y = fast_mul_mod_q2(w, T, q);
}

/////////////////// multiplying with 1/N ////////////////////////////////////

static inline void harvey_bkw_butterfly_final(uint64_t *     X, // a[j] comes from bkw_butterfly
                                              uint64_t *     Y, // a[j+N/2] comes from bkw_butterfly
                                              const mul_op_t w,
                                              const mul_op_t n_inv,
                                              const uint64_t q)
{
    const uint64_t q2 = q << 1;
    const uint64_t X1 = *X + *Y;
    const uint64_t T  = *X - *Y + q2;

    *X = fast_mul_mod_q(n_inv, X1, q);   //instead of dividing by N, we do modulus operation
    *Y = fast_mul_mod_q(w, T, q);
}


/////////function called in main//////////////////////////

void inv_ntt_ref_harvey(uint64_t       a[], // after fwd_ntt a[]
                        const uint64_t N, // N = 2^m, size of x(n)
                        const uint64_t q, // prime
                        const mul_op_t n_inv, // manual input n_inv.op found, n_inv.con not found
                        const uint64_t word_size, // 64
                        const uint64_t w[], // t->w_inv_powers.ptr
                        const uint64_t w_con[]) // t->w_inv_powers_con.ptr

{
    uint64_t t = 1;

    for(size_t m = N >> 1; m > 1; m >>= 1, t <<= 1) {
        size_t k = 0;
        for(size_t i = 0; i < m; i++) {
            const mul_op_t w1 = {w[m + i], w_con[m + i]};

            for(size_t j = k; j < k + t; j++) {
                harvey_bkw_butterfly(&a[j], &a[j + t], w1, q); // value of a[j], a[j+N/2] calculated //function defn at "fast_mul_operators.h"
            }
            k = k + (2 * t);
        }
    }

    // Final round - the harvey_bkw_butterfly, where the output is multiplies by
    // n_inv. Here m=1, k=0, t=N/2.
    const __uint128_t tmp = fast_mul_mod_q2(n_inv, w[1], q);

    // We can speed up this integer division by using barreto reduction.
    // However, as it happens only once we keep the code simple.
    const mul_op_t w1 = {tmp, (tmp << word_size) / q}; //op = manual input w or w_inverse..... con = w' = (w*2^64)/q,...for optimization, for reducing modulus step

    for(size_t j = 0; j < t; j++) {
        harvey_bkw_butterfly_final(&a[j], &a[j + t], w1, n_inv, q); // revised value of a[j], a[j+N/2] calculated,   function defn at "fast_mul_operators.h"
    }
}






