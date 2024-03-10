#include <adf.h>
#include "twiddle_last.h"

using namespace adf;

void ntt_last (input_window_int32  * restrict inputx0,
			   input_window_int32  * restrict inputx1,
			   input_window_int32  * restrict inputx2,
			   input_window_int32  * restrict inputx3,

			   output_window_int32 * restrict outputy0,
			   output_window_int32 * restrict outputy1,
			   output_window_int32 * restrict outputy2,
			   output_window_int32 * restrict outputy3){


	int32 *x0 =  (int32 *)inputx0->ptr;
	int32 *x1 =  (int32 *)inputx1->ptr;
	int32 *x2 =  (int32 *)inputx2->ptr;
	int32 *x3 =  (int32 *)inputx3->ptr;

	int32 *y0 =  (int32 *)outputy0->ptr;
	int32 *y1 =  (int32 *)outputy1->ptr;
	int32 *y2 =  (int32 *)outputy2->ptr;
	int32 *y3 =  (int32 *)outputy3->ptr;

	v8int32 * restrict tw1 = (v8int32 *) last_stage_tw1;
	v8int32 * restrict tw2 = (v8int32 *) last_stage_tw2;
	v8int32 * restrict tw3 = (v8int32 *) last_stage_tw3;

	v8int32 *xx0 = (v8int32 *)x0;
	v8int32 *xx1 = (v8int32 *)x1;
	v8int32 *xx2 = (v8int32 *)x2;
	v8int32 *xx3 = (v8int32 *)x3;

	v4int32 *yy0 = (v4int32 *)y0;
	v4int32 *yy1 = (v4int32 *)y1;
	v4int32 *yy2 = (v4int32 *)y2;
	v4int32 *yy3 = (v4int32 *)y3;

	v32int32 x1buff;

	int32 const chess_storage(%chess_alignof(int32)) tw[8] =   {1, 79, -1, -79,
											                    0,  0,   0,  0}; // twiddle for 4-point NTT (w_4)^1 == (w_4096)^1024

	v8int32 zbuff = *((v8int32 *)tw);
	v4int32 m00, m10, m20, m30, m01, m11, m21, m31;



	int32 const chess_storage(%chess_alignof(v4int32))  twt1[4] =   {1,1,1,1};

	int32 const chess_storage(%chess_alignof(v4int32))  ra[8] =   {349496, 349496, 349496, 349496, 349496, 349496, 349496, 349496}; // [( 1 << 32 ) / q]
	int32 const chess_storage(%chess_alignof(v4int32))  pa[8] =   {12289, 12289, 12289, 12289, 12289, 12289, 12289, 12289}; // q
	int32 const chess_storage(%chess_alignof(v4int32)) p12[16] =   {0, 0, 0, 0,  92253523, 92253523, 92253523, 92253523,  24578, 24578, 24578, 24578   ,92253523, 92253523, 92253523, 92253523 }; // 2q, q(1+w)

	v8int32   ra_p  = *((v8int32  * restrict)ra);
	v8int32   pa_p  = *((v8int32  * restrict)pa);
	v16int32  p12_p = *((v16int32 * restrict)p12);

	v8int32 ma0, ma1, ma2, ma3, ma1_conv, ma2_conv, ma3_conv, ma1_mod, ma2_mod, ma3_mod;
	v16int32 temp_buff, x2buff, nullv, x3buff;
	nullv = null_v16int32();
	unsigned int mask;




	for (int i = 0; i < 128; i++)
	chess_prepare_for_pipelining
	chess_loop_range(4,)
	{

		ma0 = *xx0++ ;

		ma1      = srs(mul(*xx1++, *tw1++),0);
		ma2      = srs(mul(*xx2++, *tw2++),0);
		ma3      = srs(mul(*xx3++, *tw3++),0);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);
		ma3_conv = srs(mul(srs(mul(ma3,ra_p),32),pa_p),0);

		ma1_mod  = ma1 - ma1_conv;
		ma2_mod  = ma2 - ma2_conv;
		ma3_mod  = ma3 - ma3_conv;

		x1buff = upd_w(x1buff, 0, ma0);
		x1buff = upd_w(x1buff, 1, ma1_mod);
		x1buff = upd_w(x1buff, 3, ma2_mod);
		x1buff = upd_w(x1buff, 2, ma3_mod);


		m00 = srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 0) , 0);
		m01 = srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 0) , 0);

		m10 = srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 1) , 0) ;
		m11 = srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 1) , 0) ;

		m20 = srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 2)  ,0) ;
		m21 = srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 2)  ,0) ;

		m30 = srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 3)  ,0) ;
		m31 = srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 3)  ,0) ;

		// converting negative values to positive
		temp_buff = upd_v(temp_buff, 0, m00);
		temp_buff = upd_v(temp_buff, 1, m10);
		temp_buff = upd_v(temp_buff, 2, m20);
		temp_buff = upd_v(temp_buff, 3, m30);

		mask  = lt16(temp_buff, nullv);  // checks negative element

		x2buff = add16(temp_buff , p12_p);
		x3buff = select16(mask, temp_buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98); // if x1 elem negative take elem from x2

		// modular operation
		ma1 = ext_w(x3buff, 0);
		ma2 = ext_w(x3buff, 1);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);  //m_conv = q * [(349496 * m) >> 32 ]
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;

		*yy0++ = ext_v(ma1_mod, 0);
		*yy1++ = ext_v(ma1_mod, 1);
		*yy2++ = ext_v(ma2_mod, 0);
		*yy3++ = ext_v(ma2_mod, 1);

		// for next 4 element blocks
		temp_buff = upd_v(temp_buff, 0, m01);
		temp_buff = upd_v(temp_buff, 1, m11);
		temp_buff = upd_v(temp_buff, 2, m21);
		temp_buff = upd_v(temp_buff, 3, m31);

		mask  = lt16(temp_buff, nullv);

		x2buff = add16(temp_buff , p12_p);
		x3buff = select16(mask, temp_buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98);

		ma1 = ext_w(x3buff, 0);
		ma2 = ext_w(x3buff, 1);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;

		*yy0++ = ext_v(ma1_mod, 0);
		*yy1++ = ext_v(ma1_mod, 1);
		*yy2++ = ext_v(ma2_mod, 0);
		*yy3++ = ext_v(ma2_mod, 1);

	}

}


