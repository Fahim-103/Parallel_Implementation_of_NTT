#include "ntt_last_sub.h"
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


	int32 *x0 =  (int32 *)inputx0->ptr;    //read access
	int32 *x1 =  (int32 *)inputx1->ptr;
	int32 *x2 =  (int32 *)inputx2->ptr;
	int32 *x3 =  (int32 *)inputx3->ptr;

	int32 *y0 =  (int32 *)outputy0->ptr;    //write
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

	//v4int32 *yy0 = (v4int32 *)y0;


	v8int32 m0, m1, m2, m3;
	v32int32 x1buff;

	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(int32)) tw[8] =   {1, 79, -1, -79,
											                      0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = (w_1024)^256

	v8int32 zbuff = *((v8int32 *)tw); ////////

	//modular pre-requirements
	aie::vector<int32,4> p1, p2;

	p1 = aie::broadcast<int32,4>(24578); // 2q
	p2 = aie::broadcast<int32,4>(92253523); // q(w+1) ; w is from zbuff

	v4int32 m00, m10, m20, m30, m01, m11, m21, m31;
	v4int32 mr00, mr01, mr10, mr11, mr20, mr21, mr30, mr31;


	// 8 elems in a loop; 32*8 = 256 times 4 elem ntt; 128*8 = 1024 times 4 elem ntt
	for (int i = 0; i < 128; i++){

		m0 = *xx0++ ;
		m1 = twdl_mod_mul_8(*xx1++ , *tw1++);
		m2 = twdl_mod_mul_8(*xx2++ , *tw2++);
		m3 = twdl_mod_mul_8(*xx3++ , *tw3++);

		x1buff = upd_w(x1buff, 0, m0);
		x1buff = upd_w(x1buff, 1, m1);
		x1buff = upd_w(x1buff, 3, m2);
		x1buff = upd_w(x1buff, 2, m3);


		// block of first second third fourt points, with tw multiplication

		m00 = srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 0) , 0);
		m01 = srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 0) , 0);

		m10 = srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 1) , 0) ; //first points of 4-blocks 4-point
		m11 = srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 1) , 0) ;

		m20 = srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 2)  ,0) ;
		m21 = srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 2)  ,0) ;

		m30 = srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 3)  ,0) ;
		m31 = srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 3)  ,0) ;

		//===========================================================================================================
		mr00 = twdl_mod_sub_4(m00, p1);
		mr01 = twdl_mod_sub_4(m01, p1);
		*yy0++ = mr00;
		*yy0++ = mr01;

		mr10 = twdl_mod_sub_4(m10, p2); // +p if negative and then modular
		mr11 = twdl_mod_sub_4(m11, p2);
		*yy1++ = mr10;
		*yy1++ = mr11;

		mr20 = twdl_mod_sub_4(m20, p1);
		mr21 = twdl_mod_sub_4(m21, p1);
		*yy2++ = mr20;
		*yy2++ = mr21;



		mr30 = twdl_mod_sub_4(m30, p2);
		mr31 = twdl_mod_sub_4(m31, p2);
		*yy3++ = mr30;
		*yy3++ = mr31;


	}




}


