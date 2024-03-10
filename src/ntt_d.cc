#include <adf.h>
#include "nttsub.h"
#include "include.h"
#include "twiddle.h"

using namespace adf;

void ntt_d (input_window_int32  * restrict inputx,
			output_window_int32 * restrict outputy){

    int32 *xbuff =  (int32 *)inputx->ptr;
	int32 *ybuff =  (int32 *)outputy->ptr;

	int32 * restrict tw00 = (int32 * restrict) stage0_tw;
	int32 * restrict tw01 = (int32 * restrict) stage1_tw;
	int32 * restrict tw02 = (int32 * restrict) stage2_tw;

	int32 * restrict ntt_tmp1 = ntt_1a;
	int32 * restrict ntt_tmp2 =  ntt_2a;

	stage0_ntt(xbuff, tw00, ntt_tmp1);

	for (int i = 0; i < 4; i++){

		stage1_ntt(ntt_tmp1, tw01, ntt_tmp2, i);
	}


	for (int i = 0; i < 16; i++){

		stage2_ntt(ntt_tmp2, tw02, ntt_tmp1, i);
	}

	stage3_ntt(ntt_tmp1, ybuff);

}

