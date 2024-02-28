#ifndef __NTTSUB_H__
#define __NTTSUB_H__

#include <adf.h>
#include "aie_api/aie.hpp"
#include "aie_api/utils.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/operators.hpp"


//#define rnd_floor 0

inline __attribute__((always_inline)) v4int32 twdl_mod_mul_4 ( v4int32 m, v4int32 t){   // lane by lane twidl multiplication

	aie::mask<4> msk_lt;
	aie::vector<int32,4>  mlt, mlt_con, mlt_mod, r, p;

	r=aie::broadcast<int32,4>(349496);
	p=aie::broadcast<int32,4>(12289);

	mlt = srs(mul(m, t),0);

	//modular_ntt_4 //modulus sub generation for mod = m - q [ ( m * r ) >> 32 ]
	mlt_con = srs(mul(srs(mul(mlt,r),32),p),0) ;
	mlt_mod   = mlt - mlt_con;


	return mlt_mod;
}


inline __attribute__((always_inline)) v4int32 twdl_mod_sub_4 ( v4int32 m, v4int32 pn){   // checks if the radx-4 output is negative and performs mod

	aie::mask<4> msk_lt;
	aie::vector<int32,4> m_in, m_p, msel, r, p, m_mod, m_con;

	r=aie::broadcast<int32,4>(349496);
	p=aie::broadcast<int32,4>(12289);

	// m positive or negative, i.e. m or m+p // lane selection
	m_in = m;
	m_p = m_in+ pn;
	msk_lt=aie::lt(m_in,(int32)0);
	msel=aie::select(m_in,m_p,msk_lt);


	//modular_ntt_4 //mod = m - q [ ( m * r ) >> 32 ]
	m_con = srs(mul(srs(mul(msel,r),32),p),0) ; /////////can be problem
	m_mod   = msel - m_con;

	return m_mod;
}



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

	//v4int32 * restrict tw0 = (v4int32 * restrict)  tw;
	v4int32 * restrict tw1 = (v4int32 * restrict) (tw );
	v4int32 * restrict tw2 = (v4int32 * restrict) (tw + 256);
	v4int32 * restrict tw3 = (v4int32 * restrict) (tw + 512);


	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = (w_1024)^256

	v8int32 zbuff = *((v8int32 *)tw_sub);
	//===========================================================================================================

	//modular pre-requirements
	aie::vector<int32,4> p1, p2;

	p1 = aie::broadcast<int32,4>(24578); // 2q
	p2 = aie::broadcast<int32,4>(92253523); // q(w+1) ; w is from zbuff


	v4int32 m0, m1, m2, m3;
	v4int32 mr0, mr1, mr2, mr3;
	aie::vector<int32,16> chess_storage(xa) x1buff;
	//aie::vector<int32,16> chess_storage(xd) x1buff;

	//===========================================================================================================

	x1buff = upd_v(x1buff, 0, *pi0++);
	x1buff = upd_v(x1buff, 1, *pi1++);
	x1buff = upd_v(x1buff, 3, *pi2++);
	x1buff = upd_v(x1buff, 2, *pi3++);

	//aie::print (x1buff, true, "x1buff ");

	// block of first second third fourt points, with tw multiplication

	for (int i = 0; i < 64; i++){

		m0 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; //first points of 4-blocks 4-point
		m1 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

		m2 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
		m3 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;

		//===========================================================================================================

		//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]

		mr0 = twdl_mod_sub_4(m0, p1);
		*po0++ = mr0;
		//mr0 = twdl_mod_mul_4(mr0, *t1++);  // not required for second step

		mr1 = twdl_mod_sub_4(m1, p2); // +p if negative
		mr1 = twdl_mod_mul_4(mr1, *tw1++); // tw multiplication, +1 if multp is negative
		*po1++ = mr1;

		mr2 = twdl_mod_sub_4(m2, p1);
		mr2 = twdl_mod_mul_4(mr2, *tw2++);
		*po2++ = mr2;

		mr3 = twdl_mod_sub_4(m3, p2);
		mr3 = twdl_mod_mul_4(mr3, *tw3++);
		*po3++ = mr3;


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

	//v4int32 * restrict tw0 = (v4int32 * restrict)  tw;
	v4int32 * restrict tw1 = (v4int32 * restrict) (tw );
	v4int32 * restrict tw2 = (v4int32 * restrict) (tw + 64);
	v4int32 * restrict tw3 = (v4int32 * restrict) (tw + 128);


	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = (w_1024)^256

	v8int32 zbuff = *((v8int32 *)tw_sub);
	//===========================================================================================================

	//modular pre-requirements
	aie::vector<int32,4> p1, p2;

	p1 = aie::broadcast<int32,4>(24578); // 2q
	p2 = aie::broadcast<int32,4>(92253523); // q(w+1) ; w is from zbuff


	v4int32 m0, m1, m2, m3;
	v4int32 mr0, mr1, mr2, mr3;
	aie::vector<int32,16> chess_storage(xa) x1buff;
	//aie::vector<int32,16> chess_storage(xd) x1buff;

	//===========================================================================================================

	x1buff = upd_v(x1buff, 0, *pi0++);
	x1buff = upd_v(x1buff, 1, *pi1++);
	x1buff = upd_v(x1buff, 3, *pi2++);
	x1buff = upd_v(x1buff, 2, *pi3++);

	//aie::print (x1buff, true, "x1buff ");

	// block of first second third fourt points, with tw multiplication

	for (int i = 0; i < 16; i++){

		m0 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; //first points of 4-blocks 4-point
		m1 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

		m2 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
		m3 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;

		//===========================================================================================================

		//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]

		mr0 = twdl_mod_sub_4(m0, p1);
		*po0++ = mr0;
		//mr0 = twdl_mod_mul_4(mr0, *t1++);  // not required for second step

		mr1 = twdl_mod_sub_4(m1, p2); // +p if negative
		mr1 = twdl_mod_mul_4(mr1, *tw1++); // tw multiplication, +1 if multp is negative
		*po1++ = mr1;

		mr2 = twdl_mod_sub_4(m2, p1);
		mr2 = twdl_mod_mul_4(mr2, *tw2++);
		*po2++ = mr2;

		mr3 = twdl_mod_sub_4(m3, p2);
		mr3 = twdl_mod_mul_4(mr3, *tw3++);
		*po3++ = mr3;


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

	//v4int32 * restrict tw0 = (v4int32 * restrict)  tw;
	v4int32 * restrict tw1 = (v4int32 * restrict) (tw );
	v4int32 * restrict tw2 = (v4int32 * restrict) (tw + 16);
	v4int32 * restrict tw3 = (v4int32 * restrict) (tw + 32);


	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[8] =   {1, 79, -1, -79,
											                      	  0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = (w_1024)^256

	v8int32 zbuff = *((v8int32 *)tw_sub);
	//===========================================================================================================

	//modular pre-requirements
	aie::vector<int32,4> p1, p2;

	p1 = aie::broadcast<int32,4>(24578); // 2q
	p2 = aie::broadcast<int32,4>(92253523); // q(w+1) ; w is from zbuff


	v4int32 m0, m1, m2, m3;
	v4int32 mr0, mr1, mr2, mr3;
	aie::vector<int32,16> chess_storage(xa) x1buff;
	//aie::vector<int32,16> chess_storage(xd) x1buff;

	//===========================================================================================================

	x1buff = upd_v(x1buff, 0, *pi0++);
	x1buff = upd_v(x1buff, 1, *pi1++);
	x1buff = upd_v(x1buff, 3, *pi2++);
	x1buff = upd_v(x1buff, 2, *pi3++);

	//aie::print (x1buff, true, "x1buff ");

	// block of first second third fourt points, with tw multiplication

	for (int i = 0; i < 4; i++){

		m0 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; //first points of 4-blocks 4-point
		m1 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

		m2 = srs(lmul4_sym     (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
		m3 = srs(lmul4_antisym (x1buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;

		//===========================================================================================================

		//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]

		mr0 = twdl_mod_sub_4(m0, p1);
		*po0++ = mr0;
		//mr0 = twdl_mod_mul_4(mr0, *t1++);  // not required for second step

		mr1 = twdl_mod_sub_4(m1, p2); // +p if negative
		mr1 = twdl_mod_mul_4(mr1, *tw1++);; // tw multiplication, +1 if multp is negative
		*po1++ = mr1;

		mr2 = twdl_mod_sub_4(m2, p1);
		mr2 = twdl_mod_mul_4(mr2, *tw2++);
		*po2++ = mr2;

		mr3 = twdl_mod_sub_4(m3, p2);
		mr3 = twdl_mod_mul_4(mr3, *tw3++);
		*po3++ = mr3;


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

	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw[8] =   {1, 79, -1, -79,
											                      0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = (w_1024)^256

	v8int32 zbuff = *((v8int32 *)tw); ////////

	// 4-point NTT stage00
	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[16] =   { 1, 1, 1, 1,
																		1, 7506, 4134, 7260,
																		1, 7260, 8246, 79,
																		1, 4134, 11567, 8246
																	   }; // twiddle for 16-point NTT (w_16)^1  = (w_1024)^64 // rotate for 0x2310

	v4int32 * tsub = (v4int32 *) tw_sub;


	//modular pre-requirements
	aie::vector<int32,4> p1, p2;

	p1 = aie::broadcast<int32,4>(24578); // 2q
	p2 = aie::broadcast<int32,4>(92253523); // q(w+1) ; w is from zbuff


	v4int32 m0, m1, m2, m3;
	aie::vector<int32,4> mr0, mr1, mr2, mr3;
	aie::vector<int32,16> chess_storage(xa) x1buff;
	aie::vector<int32,16> chess_storage(xd)temp_buff2;

	// update x1buff with proper permutation

/*
	x1buff = upd_v(x1buff, 0, *pi++);
	x1buff = upd_v(x1buff, 1, *pi++);
	x1buff = upd_v(x1buff, 3, *pi++);
	x1buff = upd_v(x1buff, 2, *pi++);
*/

	//aie::print(x1buff, true, "x1buff ");

	//for (int i=0; i<64;i++)
		//chess_prepare_for_pipelining
		//chess_loop_range(16,)
		//chess_unroll_loop(4)


	for(int kk = 0; kk < 4; kk++){
	        for(int jj = 0; jj < 4; jj++){
	            for(int ii = 0; ii < 4; ii++){

	            	x1buff = upd_v(x1buff, 0, *pi++);
	            	x1buff = upd_v(x1buff, 1, *pi++);
	            	x1buff = upd_v(x1buff, 3, *pi++);
	            	x1buff = upd_v(x1buff, 2, *pi++);
	            	pi += 60; // (4 + 60) = 64 vec and 64*4=256 elems

	        		//==================================================================================================

	        		m0 = srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 0) , 0);
	        		m1 = srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 1) , 0); //first points of 4-blocks 4-point
	        		m2 = srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 2)  ,0);
	        		m3 = srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 3)  ,0);

	        		//==================================================================================================
	        		//modulus sub generation for mod = t - q [ (r * t) >> 32 ]

	        		tsub++; // for mr0 all 1
	        		mr0 = twdl_mod_sub_4(m0, p1);


	        		mr1 = twdl_mod_sub_4(m1, p2); // +p if negative
	        		mr1 = twdl_mod_mul_4(mr1, *tsub++); // tw multiplication, +1 if multp is negative

	        		mr2 = twdl_mod_sub_4(m2, p1);
	        		mr2 = twdl_mod_mul_4(mr2, *tsub++);

	        		mr3 = twdl_mod_sub_4(m3, p2);
	        		mr3 = twdl_mod_mul_4(mr3, *tsub);tsub -=3;

	        		/////next stage starts/////
	        		//v16int32 x1buff, temp_buff2;

	        		x1buff = upd_v(x1buff, 0, mr0);
	        		x1buff = upd_v(x1buff, 1, mr1);
	        		x1buff = upd_v(x1buff, 2, mr2);
	        		x1buff = upd_v(x1buff, 3, mr3);

	        		temp_buff2 = aie::transpose(x1buff, 4, 4) ;

	        		//==================================================================================================
	        		// next 4-point NTT

	        		m0 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; //first points of 4-blocks 4-point
	        		m1 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) ;

	        		m2 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) ;
	        		m3 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) ;

	        		//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]

	        		mr0 = twdl_mod_sub_4(m0, p1);
	        		mr1 = twdl_mod_sub_4(m1, p2); // +p if negative
	        		mr2 = twdl_mod_sub_4(m2, p1);
	        		mr3 = twdl_mod_sub_4(m3, p2);

	        		*po = mr0[0]; po += 64;
	        		*po = mr0[1]; po += 64;
	        		*po = mr0[2]; po += 64;
	        		*po = mr0[3]; po += 64;

	        		*po = mr1[0]; po += 64;
	        		*po = mr1[1]; po += 64;
	        		*po = mr1[2]; po += 64;
	        		*po = mr1[3]; po += 64;

	        		*po = mr2[0]; po += 64;
	        		*po = mr2[1]; po += 64;
	        		*po = mr2[2]; po += 64;
	        		*po = mr2[3]; po += 64;

	        		*po = mr3[0]; po += 64;
	        		*po = mr3[1]; po += 64;
	        		*po = mr3[2]; po += 64;
	        		*po = mr3[3]; po -= 959;


	            }
	            pi -= 256; // 4-elem vectors
	            pi += 16;
	        }

	        pi -= 64;
	        pi +=4;
	    }


}


#endif
