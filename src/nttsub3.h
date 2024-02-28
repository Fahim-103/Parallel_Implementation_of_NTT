#ifndef __NTTSUB_H__
#define __NTTSUB_H__


#include <adf.h>
#include "twiddle.h"

#include "aie_api/aie.hpp"
#include "aie_api/utils.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/operators.hpp"





inline __attribute__((always_inline)) v4int32 modular_ntt ( v4int32 v){

	//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]

	int32 rr[4] = {16711935, 16711935, 16711935, 16711935};
	v4int32  r = *((v4int32*)rr);

	int32 po[8] = {257, 257, 257, 257, 257, 257, 257, 257};  // q
	v4int32 p0 = *((v4int32*)po);  // q

	v4int32 vcon, vr;
	vcon = srs(mul(srs(mul(v,r),32),p0),0) ;
	vr   = v - vcon;

	return vr;

}

inline __attribute__((always_inline)) void stage0_ntt ( int32 * restrict x, int32 * restrict y){

	v8int32 * xx = (v8int32 *) x; //16 byte 4 elem update
	v4int32 * yy = (v4int32 *) y;


	v32int32 x1buff;

	// update x1buff with proper permutation

	for (int i = 0; i < 32; i++){

		x1buff = upd_w(x1buff, 0, *xx); xx += 32;
		x1buff = upd_w(x1buff, 1, *xx); xx += 32;
		x1buff = upd_w(x1buff, 3, *xx); xx += 32;
		x1buff = upd_w(x1buff, 2, *xx); xx -= 95;
		//aie::vector<int32,32> x_vec = x1buff;
		//aie::print(x_vec, true, "x_vec ");

		// 4-point NTT twiddle

		v8int32 * zz = (v8int32 *)tw; ////////
		v8int32 zbuff = *zz;




		//modular pre-requirements
		int32 pa[4] = {514, 514, 514, 514};  // 2q
		int32 pb[4] = {4369, 4369, 4369, 4369}; // q(w+1) ; w is from zbuff

		v4int32 p1 = *((v4int32*)pa);  // 2q
		v4int32 p2 = *((v4int32*)pb);  // q(w+1) ; w is from zbuff


		//first four by four points
		// twiddl multipliers
		v4int32 * t = (v4int32 *) tww;

		v4int32 t0 = *t++;
		v4int32 t1 = *t++;
		v4int32 t2 = *t++;
		v4int32 t3 = *t;


		// block of first second third fourt points, with tw multiplication

		v4int32 m0, m1, m2, m3;

		m0 =           srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 0) , 0);
		m1 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 1) , 0) + p2 ), t1),0) ; //first points of 4-blocks 4-point
		m2 = srs(mul( (srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 2)  ,0) + p1 ), t2),0) ;
		m3 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 3)  ,0) + p2 ), t3),0) ;


		//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]
		v4int32 mr0, mr1, mr2, mr3;

		mr0 = modular_ntt(m0); *yy = mr0; yy +=32*2;
		mr1 = modular_ntt(m1); *yy = mr1; yy +=32*2;
		mr2 = modular_ntt(m2); *yy = mr2; yy +=32*2;
		mr3 = modular_ntt(m3); *yy = mr3; yy -=191;

		//next four by four points

		// twiddl multipliers
		v4int32 * t = (v4int32 *) tww;

		v4int32 t0 = *t++;
		v4int32 t1 = *t++;
		v4int32 t2 = *t++;
		v4int32 t3 = *t;


		v4int32 m00, m10, m20, m30;

		m00 =           srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 0) , 0);
		m10 = srs(mul( (srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 1) , 0) + p2 ), t1),0) ; //first points of 4-blocks 4-point
		m20 = srs(mul( (srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 2)  ,0) + p1 ), t2),0) ;
		m30 = srs(mul( (srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 3)  ,0) + p2 ), t3),0) ;


		//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]
		v4int32 mr00, mr10, mr20, mr30;

		mr00 = modular_ntt(m00); *yy = mr00; yy +=32*2;
		mr10 = modular_ntt(m10); *yy = mr10; yy +=32*2;
		mr20 = modular_ntt(m20); *yy = mr20; yy +=32*2;
		mr30 = modular_ntt(m30); *yy = mr30; yy -=191;


	}


}

inline __attribute__((always_inline)) void stage1_0_ntt ( int32 * restrict x, int32 * restrict y, int32 * restrict twd){

	v4int32 * xx = (v4int32 *) x; //16 byte 4 elem update
	v4int32 * yy = (v4int32 *) y;
	v4int32 * tt = (v4int32 *) twd;


	v16int32 x1buff;

	// update x1buff with proper permutation

	x1buff = upd_v(upd_v(x1buff, 0, *xx++), 1, *xx++);
	x1buff = upd_v(upd_v(x1buff, 3, *xx++), 2, *xx++);

	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw[8] =   {1, 256, -1, -256,
											                      0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = 16

	v8int32 zbuff = *((v8int32 *)tw); ////////

	// 4-point NTT stage00

	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[16] =   { 1,  1,      1,      1,
																		1,  64513,  16384,  65521,
																		1,  65521,  61441,  256,
																		1,  16384,  1024,   61441
																	   };

	v4int32 * t_sub = (v4int32 *) tw_sub;

	v4int32 t0_sub = *t_sub++;
	v4int32 t1_sub = *t_sub++;
	v4int32 t2_sub = *t_sub++;
	v4int32 t3_sub = *t_sub;


	//modular pre-requirements
	int32 pa[4] = {131074, 131074, 131074, 131074};  // 2q
	int32 pb[4] = {16843009, 16843009, 16843009, 16843009}; // q(w+1) ; w is from zbuff


	v4int32 p1 = *((v4int32*)pa);  // 2q
	v4int32 p2 = *((v4int32*)pb);  // q(w+1) ; w is from zbuff

	v4int32 m0, m1, m2, m3;

	m0 =           srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 0) , 0);
	m1 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 1) , 0) + p2 ), t1_sub),0) ; //first points of 4-blocks 4-point
	m2 = srs(mul( (srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 2)  ,0) + p1 ), t2_sub),0) ;
	m3 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 3)  ,0) + p2 ), t3_sub),0) ;


	//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]
	v4int32 mr0, mr1, mr2, mr3;

	mr0 = modular_ntt(m0);
	mr1 = modular_ntt(m1);
	mr2 = modular_ntt(m2);
	mr3 = modular_ntt(m3);



	/////next stage starts

	v16int32 temp_buff1, temp_buff2;

	temp_buff1 = upd_v(temp_buff1, 0, mr0);
	temp_buff1 = upd_v(temp_buff1, 1, mr1);
	temp_buff1 = upd_v(temp_buff1, 2, mr2);
	temp_buff1 = upd_v(temp_buff1, 3, mr3);

	//input vector prepare
	aie::vector<int32,16> x_vec = temp_buff1;

	//aie::print(x_vec, true, "x_vec ");

	temp_buff2 = aie::transpose(x_vec, 4, 4) ;


	// next 4-point NTT
	v4int32 k0, k1, k2, k3;

	k0 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; //first points of 4-blocks 4-point
	k1 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) + p2 ;

	k2 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) + p1 ;
	k3 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) + p2 ;


	//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]
	v4int32 kr0, kr1, kr2, kr3;

	v4int32 t0 = *tt++;
	v4int32 t1 = *tt++;
	v4int32 t2 = *tt++;
	v4int32 t3 = *tt;


	kr0 = modular_ntt(srs(mul(k0,t0),0)); //first points of 4-blocks 4-point
	kr1 = modular_ntt(srs(mul(k1,t1),0));
	kr2 = modular_ntt(srs(mul(k2,t2),0));
	kr3 = modular_ntt(srs(mul(k3,t3),0));


	*yy++ = kr0;
	*yy++ = kr1;
	*yy++ = kr2;
	*yy   = kr3;

}



inline __attribute__((always_inline)) void stage1_1_ntt ( int32 * restrict x, int32 * restrict y){

	v4int32 * xx = (v4int32 *) x; //16 byte 4 elem update
	v4int32 * yy = (v4int32 *) y;
	//v4int32 * tt = (v4int32 *) twd;


	v16int32 x1buff;

	// update x1buff with proper permutation

	x1buff = upd_v(upd_v(x1buff, 0, *xx++), 1, *xx++);
	x1buff = upd_v(upd_v(x1buff, 3, *xx++), 2, *xx++);

	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw[8] =   {1, 256, -1, -256,
											                      0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = 16

	v8int32 zbuff = *((v8int32 *)tw); ////////

	// 4-point NTT stage00

	int32 const chess_storage(%chess_alignof(v4int32)) tw_sub[16] =   { 1,  1,      1,      1,
																		1,  64513,  16384,  65521,
																		1,  65521,  61441,  256,
																		1,  16384,  1024,   61441
																	   };

	v4int32 * t_sub = (v4int32 *) tw_sub;

	v4int32 t0_sub = *t_sub++;
	v4int32 t1_sub = *t_sub++;
	v4int32 t2_sub = *t_sub++;
	v4int32 t3_sub = *t_sub;


	//modular pre-requirements
	int32 pa[4] = {131074, 131074, 131074, 131074};  // 2q
	int32 pb[4] = {16843009, 16843009, 16843009, 16843009}; // q(w+1) ; w is from zbuff


	v4int32 p1 = *((v4int32*)pa);  // 2q
	v4int32 p2 = *((v4int32*)pb);  // q(w+1) ; w is from zbuff

	v4int32 m0, m1, m2, m3;

	m0 =           srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 0) , 0);
	m1 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 1) , 0) + p2 ), t1_sub),0) ; //first points of 4-blocks 4-point
	m2 = srs(mul( (srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 2)  ,0) + p1 ), t2_sub),0) ;
	m3 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 3)  ,0) + p2 ), t3_sub),0) ;


	//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]
	v4int32 mr0, mr1, mr2, mr3;

	mr0 = modular_ntt(m0);
	mr1 = modular_ntt(m1);
	mr2 = modular_ntt(m2);
	mr3 = modular_ntt(m3);



	/////next stage starts

	v16int32 temp_buff1, temp_buff2;

	temp_buff1 = upd_v(temp_buff1, 0, mr0);
	temp_buff1 = upd_v(temp_buff1, 1, mr1);
	temp_buff1 = upd_v(temp_buff1, 2, mr2);
	temp_buff1 = upd_v(temp_buff1, 3, mr3);

	//input vector prepare
	aie::vector<int32,16> x_vec = temp_buff1;

	//aie::print(x_vec, true, "x_vec ");

	temp_buff2 = aie::transpose(x_vec, 4, 4) ;


	// next 4-point NTT
	v4int32 k0, k1, k2, k3;

	k0 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; //first points of 4-blocks 4-point
	k1 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) + p2 ;

	k2 = srs(lmul4_sym     (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) + p1 ;
	k3 = srs(lmul4_antisym (temp_buff2, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) + p2 ;


	//modulus sub generation for mod = t - q [ (16711935 * t) >> 32 ]
	v4int32 kr0, kr1, kr2, kr3;

/*
	v4int32 t0 = *tt++;
	v4int32 t1 = *tt++;
	v4int32 t2 = *tt++;
	v4int32 t3 = *tt;
*/


	kr0 = modular_ntt(k0); //first points of 4-blocks 4-point
	kr1 = modular_ntt(k1);
	kr2 = modular_ntt(k2);
	kr3 = modular_ntt(k3);


	*yy++ = kr0;
	*yy++ = kr1;
	*yy++ = kr2;
	*yy   = kr3;

}

#endif
