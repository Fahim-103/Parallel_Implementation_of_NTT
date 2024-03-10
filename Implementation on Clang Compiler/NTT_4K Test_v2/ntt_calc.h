//
// Created by secure on 11/7/23.
//

#ifndef NTT_CALC_H
#define NTT_CALC_H


#include <bits/stdint-uintn.h>
#include <glob.h>
#include <stdlib.h>
#include <stdio.h>

///////////////////////begining///////////////////////


#define WORD_SIZE 64 //no impact when parameters are not long integers

typedef struct mul_op_s {
    __uint128_t op;  //w =
    __uint128_t con; //w converted///////////////////////////////////////////////
} mul_op_t;

//structure input

typedef struct aligned64_ptr_s {
    void *    base;
    uint64_t *ptr;
} aligned64_ptr_t;

//structure test case creation

typedef struct test_case_s {
    // These parameters are predefined
    uint64_t m; // n= 2^m
    uint64_t q; //ring, or, field size Rq, q = lN + 1
    uint64_t w; //primitive Nth root OF unity
    uint64_t w_inv; // manual input w^(-1) mod q
    mul_op_t n_inv; //

    // These parameters are dinamically computed based on the above values.
    uint64_t        n; // n= 2^m
    uint64_t        qneg;
    uint64_t        q2; // q2 = 2q
    uint64_t        q4; // q4 = 4q


    uint64_t   w_powers[4096]; // an array w^(nk), allocate_aligned_array(), calc_w()
    uint64_t   w_powers_con[4096]; //an array w'^(nk)
    uint64_t   w_inv_powers[4096]; //this is required,, calc_w_inv(),
    uint64_t   w_inv_powers_con[4096];


} test_case_t;


static inline uint64_t bit_rev_idx(uint64_t idx, uint64_t width)
{
    uint64_t ret = 0;
    while(width > 0) {
        width--;
        ret += ((idx & 1) << width);
        idx >>= 1;
    }

    return ret;
}

static inline void bit_rev(uint64_t       w_powers[],
                           const uint64_t w[],
                           const uint64_t N,
                           const uint64_t width)
{
    for(size_t i = 0; i < N; i++) {
        w_powers[i] = w[i];
    }
}

static inline void calc_w(uint64_t       w_powers_rev[],  // value goes t->w_powers.ptr, w_powers.ptr = somthing + base
                          const uint64_t w, // nth rootOFunity manual input
                          const uint64_t N, // x(n) size, n = 2^m
                          const uint64_t q, // prime number
                          const uint64_t width) //calculating N numbers w_powers

{
    uint64_t w_powers[N];
    w_powers[0] = 1;
    for(size_t i = 1; i < N; i++) {
        w_powers[i] = (uint64_t)(((__uint128_t)w_powers[i - 1] * w) % q); //wpowers[i] = [w^(nk)] % q
        //printf("w_powers[%ld] = %ld\n", i, w_powers[i]);
    }

    bit_rev(w_powers_rev, w_powers, N, width); // clclated w_powers[] are bit_rev, set at correspond



}

static inline void calc_w_inv(uint64_t       w_inv_rev[],
                              const uint64_t w_inv,
                              const uint64_t N,
                              const uint64_t q,
                              const uint64_t width)
{
    uint64_t w_inv_powers[N];
    w_inv_powers[0] = 1;
    for(size_t i = 1; i < N; i++) {
        w_inv_powers[i] = (uint64_t)(((__uint128_t)w_inv_powers[i - 1] * w_inv) % q);
    }

    bit_rev(w_inv_rev, w_inv_powers, N, width);
}
static inline void calc_w_con(uint64_t       w_con[],  //w'
                              const uint64_t w[],  //w
                              const uint64_t N, //2^m
                              const uint64_t q, // prime
                              const uint64_t word_size) //beta, B = 2^32
{
    for(size_t i = 0; i < N; i++) {
        w_con[i] = ((__uint128_t)w[i] << word_size) / q;   //w' = [(w*(2^32))/q]
        //printf("calc w_p %lu w_con %lu\n", w[i], w_con[i]);
    }
}

static uint64_t
calc_ninv_con(const uint64_t Ninv, const uint64_t q, const uint64_t word_size)
{
    return ((__uint128_t)Ninv << word_size) / q;  //(Ninv*(2^64)/q
}

///////////////////////function required for calc_w,.........closed.........


static inline int _init_test(test_case_t *t) // called at main.c - _init_tes_init_test_cases()
// takes manual_input of m, q, w, w_inv, n_inv.op
{
    // For brevity
    const uint64_t q     = t->q;
    const uint64_t w     = t->w;
    const uint64_t m     = t->m;
    const uint64_t w_inv = t->w_inv;

    ///////// calculate n//////
    const uint64_t n     = 1UL << t->m; //  rshft to m bits, 2^m

    //printf("nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn %lu\n",n);

    t->n         = n;


    // Prepare radix-2 w-powers
    //allocate_aligned_array(&t->w_powers, n);
    calc_w(t->w_powers, w, n, q, m);  //check the destination w_powers.ptr
                                                            // returns array of w^powers after bit reversed order

    //allocate_aligned_array(&t->w_powers_con, n);
    calc_w_con(t->w_powers_con, t->w_powers, n, q, WORD_SIZE); //check the destination w_powers_con.ptr

   // allocate_aligned_array(&t->w_inv_powers, n);
    calc_w_inv(t->w_inv_powers, w_inv, n, q, m);

   // allocate_aligned_array(&t->w_inv_powers_con, n);
    calc_w_con(t->w_inv_powers_con, t->w_inv_powers, n, q, WORD_SIZE);

    t->n_inv.con = calc_ninv_con(t->n_inv.op, q, WORD_SIZE);
    t->q2        = 2 * q;
    t->q4        = 4 * q;


    return 1;
}


///////////////////////CLOSED test_case input vector and corresponding paramaeter preparation CLOSED///////////




#endif //UNTITLED1_NTT_CALC_H
