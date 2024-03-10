#include <stdio.h>
//#include "ntt_calc.h"
#include "ntt_fwd_inv.h"

int main() {

    //static test_case_t tests = {.m = 12, .q = 0x10001, .w = 21, .w_inv = 49933, .n_inv.op = 65521};

   static test_case_t tests = {.m = 12, .q = 12289, .w = 11, .w_inv = 5586, .n_inv.op = 12286};
   // static test_case_t tests = {.m = 10, .q = 12289, .w = 11, .w_inv = 5586, .n_inv.op = 12277};


    //static test_case_t tests = {.m = 4, .q = 257, .w = 32, .w_inv = 2305, .n_inv.op = 755};

    _init_test(&tests);  //test_case parameter initiation

    const test_case_t *t;
    t = &tests;


    /**************************************************************************/

    uint64_t a_in[t->n]; uint64_t * a_i = a_in;
    uint64_t a[t->n];

    for (int i = 0; i<(t->n); i++){

        a_in[i] = i;
        a[i] = i;
    }

    printf("Printing input\n");

    for (int i = 0; i < t->n; i++) {
        printf("%ld ", a_in[i]);
    }
    printf("\n");


    //generate Forward NTT -- a

    fwd_ntt_ref_harvey(a, t->n, t->q, t->w_powers, t->w_powers_con, t->m);

/*

    uint64_t a_fwd[t->n];
    for(int i=0; i< t->n; i++){

        uint64_t tt = bit_rev_idx(i,t->m);
     //   printf("%d to %d\n", i, tt);
        a_fwd[tt] = a[i];

    }

*/


    printf("Printing Forward NTT of a bit_rev\n");

    for (int i = 0; i < t->n; i++) {
        printf("%ld ", a[i]);
    }

    printf("\n");



    //generate Inverse NTT

    inv_ntt_ref_harvey(a, t->n, t->q, t->n_inv, WORD_SIZE, t->w_inv_powers,
                       t->w_inv_powers_con, t->m);


    //////////////////////////////////////////////////////////

    printf("Printing Reverse NTT of a bit_rev\n");

    for (int i = 0; i < t->n; i++) {
        printf("%ld ", a[i]);
    }

    printf("\n");



    //////////////////////////////////            n1*n2              /////////////////////////////////////////////////

    int n1 = 4;  //////
    int m1 = 2;  //////

    int n2 = (t->n)/n1;  //////
    int m2 = t->m - m1;  ///////

    //////perform n2-point NTT creating arrays/////

    uint64_t a1[n1];

    uint64_t a2_0[n2]; uint64_t * p0 = a2_0;
    uint64_t a2_1[n2]; uint64_t * p1 = a2_1;
    uint64_t a2_2[n2]; uint64_t * p2 = a2_2;
    uint64_t a2_3[n2]; uint64_t * p3 = a2_3;

    for (int i =0; i< (t->n)/n1; i++){

        *p0++ = *a_i++;
        *p1++ = *a_i++;
        *p2++ = *a_i++;
        *p3++ = *a_i++;
    }

    printf("1st vector: ");

    for (int i =0; i< (t->n)/n1; i++){

        printf( "%lu ", a2_0[i]);
    }

    printf("\n");

    //////////////////////////////////////////////////////////////

    printf("2nd vector: ");
    for (int i =0; i< (t->n)/n1; i++){

        printf("%lu ", a2_1[i]);
    }

    printf("\n");

    //////////////////////////////////////////////////////////////

    printf("3rd vector: ");
    for (int i =0; i< (t->n)/n1; i++){

        printf("%lu ", a2_2[i]);
    }

    printf("\n");

    //////////////////////////////////////////////////////////////

    printf("4th vector: ");
    for (int i =0; i< (t->n)/n1; i++){

        printf("%lu ", a2_3[i]);
    }

    printf("\n");

    //////////////////////////////////////////////////////////////

    uint64_t w_powers_n2[n2];
    uint64_t w_powers_con_n2[n2];

    for(int i = 0; i< n2; i++){

        w_powers_n2[i]  = t->w_powers[i*n1];
        w_powers_con_n2[i] = t->w_powers_con[i*n1];
    }

    fwd_ntt_ref_harvey(a2_0, n2, t->q, w_powers_n2, w_powers_con_n2, m2);
    fwd_ntt_ref_harvey(a2_1, n2, t->q, w_powers_n2, w_powers_con_n2, m2);
    fwd_ntt_ref_harvey(a2_2, n2, t->q, w_powers_n2, w_powers_con_n2, m2);
    fwd_ntt_ref_harvey(a2_3, n2, t->q, w_powers_n2, w_powers_con_n2, m2);


    printf("3rd fwd vector: ");
    for (int i =0; i< (t->n)/n1; i++){

        printf("%lu ", a2_2[i]);
    }

    printf("\n");



    /////////twiddle multiplication////////

    int j =0;
        for (int i = 0; i < (t->n) / n1; i++) {
            a2_0[i] = (a2_0[i] * t->w_powers[i * j]) % t->q;  // i means column
        }
    j++;
        for (int i = 0; i < (t->n) / n1; i++) {
            a2_1[i] = (a2_1[i] * t->w_powers[i * j]) % t->q;  // i means column
        }
    j++;
    for (int i = 0; i < (t->n) / n1; i++) {
        a2_2[i] = (a2_2[i] * t->w_powers[i * j]) % t->q;  // i means column
    }
    j++;
    for (int i = 0; i < (t->n) / n1; i++) {
        a2_3[i] = (a2_3[i] * t->w_powers[i * j]) % t->q;  // i means column
    }


    //////////////////////////////////////////////

    uint64_t a3[n1]; //uint64_t * pp = a3;


    for(int i = 0; i < n2; i++){

        a3[0] = a2_0[i];
        a3[1] = a2_1[i];
        a3[2] = a2_2[i];
        a3[3] = a2_3[i];

        uint64_t w_powers_n1[n1];
        uint64_t w_powers_con_n1[n1];


        for(int i = 0; i< n1; i++){

            w_powers_n1[i]  = t->w_powers[i*n2];
            w_powers_con_n1[i] = t->w_powers_con[i*n2];
        }

        fwd_ntt_ref_harvey(a3, n1, t->q, w_powers_n1, w_powers_con_n1, m1);


        /*printf("test: ");
        for (int i =0; i< n1; i++){

            printf("%lu ", a3[i]);
        }

        printf("\n");*/




    }





    return 0;
}
