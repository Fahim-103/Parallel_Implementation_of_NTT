#ifndef __NTTSUB_H__
#define __NTTSUB_H__

inline __attribute__((always_inline)) void stage0_ntt (int32 * restrict x, int32 * restrict tw,
														  int32 * restrict y){


	v4int32 * restrict pi0 = (v4int32 * restrict)  x;
	v4int32 * restrict pi1 = (v4int32 * restrict) (x + 256);
	v4int32 * restrict pi2 = (v4int32 * restrict) (x + 512);
	v4int32 * restrict pi3 = (v4int32 * restrict) (x + 768);


	v4int32 * restrict po0 = (v4int32 * restrict)  y;
	v4int32 * restrict po1 = (v4int32 * restrict) (y + 256);
	v4int32 * restrict po2 = (v4int32 * restrict) (y + 512);
	v4int32 * restrict po3 = (v4int32 * restrict) (y + 768);


	v4int32 * restrict tw1 = (v4int32 * restrict) (tw );
	v4int32 * restrict tw2 = (v4int32 * restrict) (tw + 256);
	v4int32 * restrict tw3 = (v4int32 * restrict) (tw + 512);


	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0}; //zbuff twiddle for 4-point NTT (w_4)^1 == (w_4096)^1024

	v8int32 zbuff = *((v8int32 * restrict)tw_sub);
	v4int32 m0, m1, m2, m3;


	// pre-requirements for modular operation

	int32 const chess_storage(%chess_alignof(v4int32)) twt1[4] =   {1,1,1,1};
	int32 const chess_storage(%chess_alignof(v4int32))   ra[8] =   {349496, 349496, 349496, 349496, 349496, 349496, 349496, 349496}; // [( 1 << 32 ) / q]
	int32 const chess_storage(%chess_alignof(v4int32))   pa[8] =   {12289, 12289, 12289, 12289, 12289, 12289, 12289, 12289}; // q
	int32 const chess_storage(%chess_alignof(v4int32)) p12[16] =   {0, 0, 0, 0,  92253523, 92253523, 92253523, 92253523,  24578, 24578, 24578, 24578   ,92253523, 92253523, 92253523, 92253523 }; // 2q, q(1+w)

	v4int32   twt1a = *((v4int32  * restrict)twt1);
	v8int32   ra_p  = *((v8int32  * restrict)ra);
	v8int32   pa_p  = *((v8int32  * restrict)pa);
	v16int32  p12_p = *((v16int32 * restrict)p12);

	v8int32 ma1, ma2, ma1_conv, ma2_conv, ma1_mod, ma2_mod, twa1, twa2;
	v16int32 x1buff, x2buff, nullv, x3buff;
	nullv = null_v16int32();
	unsigned int mask;


	x1buff = upd_v(x1buff, 0, *pi0++);
	x1buff = upd_v(x1buff, 1, *pi1++);
	x1buff = upd_v(x1buff, 3, *pi2++);
	x1buff = upd_v(x1buff, 2, *pi3++);

	for (int i = 0; i < 64; i++)
	chess_prepare_for_pipelining
	chess_loop_range(4,)
	{

		m0 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; // four first-points of the four-blocks
		m1 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ; // four secnd-points of the four-blocks

		m2 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
		m3 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;


		// converting negative values to positive

		x1buff = upd_v(x1buff, 0, m0);
		x1buff = upd_v(x1buff, 1, m1);
		x1buff = upd_v(x1buff, 2, m2);
		x1buff = upd_v(x1buff, 3, m3);

		mask  = lt16(x1buff, nullv);  // checks negative element

		x2buff = add16(x1buff , p12_p);
		x3buff = select16(mask, x1buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98); // if x1 elem negative, take elem from x2

		// modular operation

		ma1 = ext_w(x3buff, 0);
		ma2 = ext_w(x3buff, 1);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);  //m_conv = q * [(349496 * m) >> 32 ]
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;


		// twiddle multiplication & afterwards modular operation

		twa1 =  upd_v(twa1, 0, twt1a);
		twa1 =  upd_v(twa1, 1, *tw1++);
		twa2 =  upd_v(twa2, 0, *tw2++);
		twa2 =  upd_v(twa2, 1, *tw3++);


		ma1 = srs(mul(ma1_mod, twa1),0);
		ma2 = srs(mul(ma2_mod, twa2),0);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;

		*po0++ = ext_v(ma1_mod, 0);
		*po1++ = ext_v(ma1_mod, 1);
		*po2++ = ext_v(ma2_mod, 0);
		*po3++ = ext_v(ma2_mod, 1);


		x1buff = upd_v(x1buff, 0, *pi0++);
		x1buff = upd_v(x1buff, 1, *pi1++);
		x1buff = upd_v(x1buff, 3, *pi2++);
		x1buff = upd_v(x1buff, 2, *pi3++);

	}

}


inline __attribute__((always_inline)) void stage1_ntt (int32 * restrict x, int32 * restrict tw,
														  int32 * restrict y, int i){


	v4int32 * restrict pi0 = (v4int32 * restrict)  (x + i*256);
	v4int32 * restrict pi1 = (v4int32 * restrict) ((x + i*256) + 64);
	v4int32 * restrict pi2 = (v4int32 * restrict) ((x + i*256) + 128);
	v4int32 * restrict pi3 = (v4int32 * restrict) ((x + i*256) + 192);


	v4int32 * restrict po0 = (v4int32 * restrict)  (y + i*256);
	v4int32 * restrict po1 = (v4int32 * restrict) ((y + i*256) + 64);
	v4int32 * restrict po2 = (v4int32 * restrict) ((y + i*256) + 128);
	v4int32 * restrict po3 = (v4int32 * restrict) ((y + i*256) + 192);

	v4int32 * restrict tw1 = (v4int32 * restrict) (tw );
	v4int32 * restrict tw2 = (v4int32 * restrict) (tw + 64);
	v4int32 * restrict tw3 = (v4int32 * restrict) (tw + 128);

	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0}; //zbuff twiddle for 4-point NTT (w_4)^1 == (w_4096)^1024

	v8int32 zbuff = *((v8int32 * restrict)tw_sub);
	v4int32 m0, m1, m2, m3;

	int32 const chess_storage(%chess_alignof(v4int32))  twt1[4] =   {1,1,1,1};
	int32 const chess_storage(%chess_alignof(v4int32))  ra[8] =   {349496, 349496, 349496, 349496, 349496, 349496, 349496, 349496};
	int32 const chess_storage(%chess_alignof(v4int32))  pa[8] =   {12289, 12289, 12289, 12289, 12289, 12289, 12289, 12289};
	int32 const chess_storage(%chess_alignof(v4int32)) p12[16] =   {0, 0, 0, 0,  92253523, 92253523, 92253523, 92253523,  24578, 24578, 24578, 24578   ,92253523, 92253523, 92253523, 92253523 };

	v4int32   twt1a = *((v4int32  * restrict)twt1);
	v8int32   ra_p  = *((v8int32  * restrict)ra);
	v8int32   pa_p  = *((v8int32  * restrict)pa);
	v16int32  p12_p = *((v16int32 * restrict)p12);

	v8int32 ma1, ma2, ma1_conv, ma2_conv, ma1_mod, ma2_mod, twa1, twa2;
	v16int32 x1buff, x2buff, nullv, x3buff;
	nullv = null_v16int32();
	unsigned int mask;

	x1buff = upd_v(x1buff, 0, *pi0++);
	x1buff = upd_v(x1buff, 1, *pi1++);
	x1buff = upd_v(x1buff, 3, *pi2++);
	x1buff = upd_v(x1buff, 2, *pi3++);

	for (int i = 0; i < 16; i++)
	chess_prepare_for_pipelining
	chess_loop_range(4,)
	{

		m0 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ;
		m1 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

		m2 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
		m3 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;


		x1buff = upd_v(x1buff, 0, m0);
		x1buff = upd_v(x1buff, 1, m1);
		x1buff = upd_v(x1buff, 2, m2);
		x1buff = upd_v(x1buff, 3, m3);

		mask  = lt16(x1buff, nullv);
		x2buff = add16(x1buff , p12_p);
		x3buff = select16(mask, x1buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98);

		ma1 = ext_w(x3buff, 0);
		ma2 = ext_w(x3buff, 1);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;


		twa1 =  upd_v(twa1, 0, twt1a);
		twa1 =  upd_v(twa1, 1, *tw1++);
		twa2 =  upd_v(twa2, 0, *tw2++);
		twa2 =  upd_v(twa2, 1, *tw3++);

		ma1 = srs(mul(ma1_mod, twa1),0);
		ma2 = srs(mul(ma2_mod, twa2),0);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;

		*po0++ = ext_v(ma1_mod, 0);
		*po1++ = ext_v(ma1_mod, 1);
		*po2++ = ext_v(ma2_mod, 0);
		*po3++ = ext_v(ma2_mod, 1);

		x1buff = upd_v(x1buff, 0, *pi0++);
		x1buff = upd_v(x1buff, 1, *pi1++);
		x1buff = upd_v(x1buff, 3, *pi2++);
		x1buff = upd_v(x1buff, 2, *pi3++);


	}

}

inline __attribute__((always_inline)) void stage2_ntt (int32 * restrict x, int32 * restrict tw,
														  int32 * restrict y, int i){


	v4int32 * restrict pi0 = (v4int32 * restrict)  (x + i*64);
	v4int32 * restrict pi1 = (v4int32 * restrict) ((x + i*64) + 16);
	v4int32 * restrict pi2 = (v4int32 * restrict) ((x + i*64) + 32);
	v4int32 * restrict pi3 = (v4int32 * restrict) ((x + i*64) + 48);


	v4int32 * restrict po0 = (v4int32 * restrict)  (y + i*64);
	v4int32 * restrict po1 = (v4int32 * restrict) ((y + i*64) + 16);
	v4int32 * restrict po2 = (v4int32 * restrict) ((y + i*64) + 32);
	v4int32 * restrict po3 = (v4int32 * restrict) ((y + i*64) + 48);

	v4int32 * restrict tw1 = (v4int32 * restrict) (tw );
	v4int32 * restrict tw2 = (v4int32 * restrict) (tw + 16);
	v4int32 * restrict tw3 = (v4int32 * restrict) (tw + 32);


	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0}; //zbuff twiddle for 4-point NTT (w_4)^1 == (w_4096)^1024


	v8int32 zbuff = *((v8int32 * restrict)tw_sub);
	v4int32 m0, m1, m2, m3;

	int32 const chess_storage(%chess_alignof(v4int32))  twt1[4] =   {1,1,1,1};
	int32 const chess_storage(%chess_alignof(v4int32))  ra[8]   =   {349496, 349496, 349496, 349496, 349496, 349496, 349496, 349496};
	int32 const chess_storage(%chess_alignof(v4int32))  pa[8]   =   {12289, 12289, 12289, 12289, 12289, 12289, 12289, 12289};
	int32 const chess_storage(%chess_alignof(v4int32))  p12[16] =   {0, 0, 0, 0,  92253523, 92253523, 92253523, 92253523,  24578, 24578, 24578, 24578   ,92253523, 92253523, 92253523, 92253523 };

	v8int32   ra_p  = *((v8int32  * restrict)ra);
	v8int32   pa_p  = *((v8int32  * restrict)pa);
	v4int32   twt1a = *((v4int32  * restrict)twt1);
	v16int32  p12_p = *((v16int32 * restrict)p12);

	v8int32 ma1, ma2, ma1_conv, ma2_conv, ma1_mod, ma2_mod, twa1, twa2;
	v16int32 x1buff, x2buff, nullv, x3buff;
	nullv = null_v16int32();
	unsigned int mask;

	x1buff = upd_v(x1buff, 0, *pi0++);
	x1buff = upd_v(x1buff, 1, *pi1++);
	x1buff = upd_v(x1buff, 3, *pi2++);
	x1buff = upd_v(x1buff, 2, *pi3++);


	for (int i = 0; i < 4; i++)
	chess_prepare_for_pipelining
	chess_loop_range(4,)
	{

		m0 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ;
		m1 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

		m2 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
		m3 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;


		x1buff = upd_v(x1buff, 0, m0);
		x1buff = upd_v(x1buff, 1, m1);
		x1buff = upd_v(x1buff, 2, m2);
		x1buff = upd_v(x1buff, 3, m3);

		mask  = lt16(x1buff, nullv);

		x2buff = add16(x1buff , p12_p);
		x3buff = select16(mask, x1buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98);

		ma1 = ext_w(x3buff, 0);
		ma2 = ext_w(x3buff, 1);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;

		twa1 =  upd_v(twa1, 0, twt1a);
		twa1 =  upd_v(twa1, 1, *tw1++);
		twa2 =  upd_v(twa2, 0, *tw2++);
		twa2 =  upd_v(twa2, 1, *tw3++);


		ma1 = srs(mul(ma1_mod, twa1),0);
		ma2 = srs(mul(ma2_mod, twa2),0);

		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

		ma1_mod = ma1 - ma1_conv;
		ma2_mod = ma2 - ma2_conv;

		*po0++ = ext_v(ma1_mod, 0);
		*po1++ = ext_v(ma1_mod, 1);
		*po2++ = ext_v(ma2_mod, 0);
		*po3++ = ext_v(ma2_mod, 1);

		x1buff = upd_v(x1buff, 0, *pi0++);
		x1buff = upd_v(x1buff, 1, *pi1++);
		x1buff = upd_v(x1buff, 3, *pi2++);
		x1buff = upd_v(x1buff, 2, *pi3++);

	}

}

inline __attribute__((always_inline)) void stage3_ntt (int32 * restrict x,
													   int32 * restrict y){


	v4int32 * restrict pi = (v4int32 * restrict) x;
	int32 * restrict po = (int32 * restrict) y;

	int32 const chess_storage(%chess_alignof(v4int32)) tw[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0};

	v8int32 zbuff = *((v8int32 *)tw);

	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[16] =   { 1, 1, 1, 1,
																		1, 7506, 4134, 7260,
																		1, 7260, 8246, 79,
																		1, 4134, 11567, 8246
																	   }; // twiddle for 16-point NTT (w_16)^1  == (w_4096)^256 // rotate for 0x2310

	v4int32 * tsub = (v4int32 * restrict) tw_sub;
	v4int32 m0, m1, m2, m3;


	int32 const chess_storage(%chess_alignof(v4int32))  twt1[4] =   {1,1,1,1};
	int32 const chess_storage(%chess_alignof(v4int32))  ra[8] =   {349496, 349496, 349496, 349496, 349496, 349496, 349496, 349496};
	int32 const chess_storage(%chess_alignof(v4int32))  pa[8] =   {12289, 12289, 12289, 12289, 12289, 12289, 12289, 12289};
	int32 const chess_storage(%chess_alignof(v4int32))  p12[16] =   {0, 0, 0, 0,  92253523, 92253523, 92253523, 92253523,  24578, 24578, 24578, 24578   ,92253523, 92253523, 92253523, 92253523 };


	v4int32   twt1a = *((v4int32  * restrict)twt1);
	v8int32   ra_p  = *((v8int32  * restrict)ra);
	v8int32   pa_p  = *((v8int32  * restrict)pa);
	v16int32  p12_p = *((v16int32 * restrict)p12);

	v8int32 ma1, ma2, ma1_conv, ma2_conv, ma1_mod, ma2_mod, twa1, twa2;
	v16int32 x1buff, x2buff, nullv, x3buff, temp_buff2;
	nullv = null_v16int32();
	unsigned int mask;

	for(int kk = 0; kk < 4; kk++){
	        for(int jj = 0; jj < 4; jj++){
	            for(int ii = 0; ii < 4; ii++)
	            chess_prepare_for_pipelining
	            chess_loop_range(4,)
	            {

	            	x1buff = upd_v(x1buff, 0, *pi++);
	            	x1buff = upd_v(x1buff, 1, *pi++);
	            	x1buff = upd_v(x1buff, 3, *pi++);
	            	x1buff = upd_v(x1buff, 2, *pi++);
	            	pi += 60; // (4 + 60) = 64 vec and 64*4=256 elems

	        		m0 = srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 0) , 0);
	        		m1 = srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 1) , 0);
	        		m2 = srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 2)  ,0);
	        		m3 = srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 3)  ,0);

	        		x1buff = upd_v(x1buff, 0, m0);
	        		x1buff = upd_v(x1buff, 1, m1);
	        		x1buff = upd_v(x1buff, 2, m2);
	        		x1buff = upd_v(x1buff, 3, m3);

	        		mask  = lt16(x1buff, nullv);

	        		x2buff = add16(x1buff , p12_p);
	        		x3buff = select16(mask, x1buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98);

	        		ma1 = ext_w(x3buff, 0);
	        		ma2 = ext_w(x3buff, 1);

	        		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
	        		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

	        		ma1_mod = ma1 - ma1_conv;
	        		ma2_mod = ma2 - ma2_conv;

	        		twa1 =  upd_v(twa1, 0, *tsub);
	        		twa1 =  upd_v(twa1, 1, *(tsub + 1));
	        		twa2 =  upd_v(twa2, 0, *(tsub + 2));
	        		twa2 =  upd_v(twa2, 1, *(tsub + 3));


	        		ma1 = srs(mul(ma1_mod, twa1),0);
	        		ma2 = srs(mul(ma2_mod, twa2),0);

	        		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
	        		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

	        		ma1_mod = ma1 - ma1_conv;
	        		ma2_mod = ma2 - ma2_conv;

	        		// next stage starts

	        		x1buff = upd_w(x1buff, 0, ma1_mod);
	        		x1buff = upd_w(x1buff, 1, ma2_mod);

	        		temp_buff2 = shuffle16(x1buff,0,0xD951C840,0xFB73EA62);  // transpose

	        		m0 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ;
	        		m1 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

	        		m2 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
	        		m3 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;

	        		x1buff = upd_v(x1buff, 0, m0);
	        		x1buff = upd_v(x1buff, 1, m1);
	        		x1buff = upd_v(x1buff, 2, m2);
	        		x1buff = upd_v(x1buff, 3, m3);

	        		mask  = lt16(x1buff, nullv);

	        		x2buff = add16(x1buff , p12_p);
	        		x3buff = select16(mask, x1buff, 0, 0x76543210, 0xFEDCBA98, x2buff , 0, 0x76543210, 0xFEDCBA98);

	        		ma1 = ext_w(x3buff, 0);
	        		ma2 = ext_w(x3buff, 1);

	        		ma1_conv = srs(mul(srs(mul(ma1,ra_p),32),pa_p),0);
	        		ma2_conv = srs(mul(srs(mul(ma2,ra_p),32),pa_p),0);

	        		ma1_mod = ma1 - ma1_conv;
	        		ma2_mod = ma2 - ma2_conv;

	        		*po = ext_elem (ma1_mod, 0); po += 64;
	        		*po = ext_elem (ma1_mod, 1); po += 64;
	        		*po = ext_elem (ma1_mod, 2); po += 64;
	        		*po = ext_elem (ma1_mod, 3); po += 64;

	        		*po = ext_elem (ma1_mod, 4); po += 64;
	        		*po = ext_elem (ma1_mod, 5); po += 64;
	        		*po = ext_elem (ma1_mod, 6); po += 64;
	        		*po = ext_elem (ma1_mod, 7); po += 64;

	        		*po = ext_elem (ma2_mod, 0); po += 64;
	        		*po = ext_elem (ma2_mod, 1); po += 64;
	        		*po = ext_elem (ma2_mod, 2); po += 64;
	        		*po = ext_elem (ma2_mod, 3); po += 64;

	        		*po = ext_elem (ma2_mod, 4); po += 64;
	        		*po = ext_elem (ma2_mod, 5); po += 64;
	        		*po = ext_elem (ma2_mod, 6); po += 64;
	        		*po = ext_elem (ma2_mod, 7); po -= 959;


	            }
	            pi -= 256;
	            pi += 16;
	        }

	        pi -= 64;
	        pi +=4;
	    }

}




#endif
