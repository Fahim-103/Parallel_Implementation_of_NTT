#ifndef __NTTSUB_H__
#define __NTTSUB_H__


#include <adf.h>
#include "twiddle.h"

#include "aie_api/aie.hpp"
#include "aie_api/utils.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/operators.hpp"


inline __attribute__((always_inline)) void stage0_ntt ( int32 * restrict x, int32 * restrict y){

	v4int32 * xx = (v4int32 *) x; //16 byte 4 elem update
	v8int32 * yy = (v8int32 *) y;
	v16int32 x1buff, x2buff;

	// update x1buff with proper permutation

	x1buff = upd_v(upd_v(x1buff, 0, *xx++), 1, *xx++);
	x1buff = upd_v(upd_v(x1buff, 3, *xx++), 2, *xx++);

	// 4-point NTT twiddle

	v8int32 * zz = (v8int32 *)tw; ////////
	v8int32 zbuff = *zz;

	// 4-point NTT stage00

	v4int32 m0, m1, m2, m3;

	v4int32 * t = (v4int32 *) tww;

	v4int32 t0 = *t++;
	v4int32 t1 = *t++;
	v4int32 t2 = *t++;
	v4int32 t3 = *t;


	int32 pa[4] = {514, 514, 514, 514};  // 2q
	int32 pb[4] = {4369, 4369, 4369, 4369}; // q(w+1) ; w is from zbuff


	v4int32 p1 = *((v4int32*)pa);  // 2q
	v4int32 p2 = *((v4int32*)pb);  // q(w+1) ; w is from zbuff



	// n0 = srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 1) , 0);srs(mul ( (n0+p2) , t1) , 0)

	m0 = srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 0) , 0);  x2buff = upd_v(x2buff, 0, m0);
	m1 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 1) , 0) + p2 ), t1),0) ; x2buff = upd_v(x2buff, 1, m1);
	m2 = srs(mul( (srs(lmul4_sym     (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 2)  ,0) + p1 ), t2),0) ; x2buff = upd_v(x2buff, 2, m2);
	m3 = srs(mul( (srs(lmul4_antisym (x1buff, 0, 0x2310, 4, 12, zbuff, 0, 0x0000, 3)  ,0) + p2 ), t3),0) ; x2buff = upd_v(x2buff, 3, m3);


	// modulus ops
	/// only uint8 accepted for vector



	int32 po[8] = {257, 257, 257, 257, 257, 257, 257, 257};  // q
	v8int32 p0 = *((v8int32*)po);  // q


	aie::vector<int32,8> x_vec_first =   ext_w(x2buff,0);  // '0' means first 8 lanes; w means 256 bits
	aie::vector<int32,8> x_vec_second =  ext_w(x2buff,1);


	///modulus sub generation for mod = t - q [ (255 * t) >> 16 ]
	aie::vector<int32,8> x_con_f = (srs(aie::mul(x_vec_first, 255),0)); // mul saved in acc80, v16acc80 not possible; no srs() for v16int32
	v8int32 x_con_f_down_q = srs( mul(aie::downshift(x_con_f,16), p0), 0); // bit and 0 for first bit unsigned

	aie::vector<int32,8> x_con_s = (srs(aie::mul(x_vec_second, 255),0));
	v8int32 x_con_s_down_q = srs( mul(aie::downshift(x_con_s,16), p0), 0);

	///modulus sub execution
	v16int32 x2nbuff;
	x2nbuff= upd_w(x2nbuff, 0, x_con_f_down_q);
	x2nbuff = upd_w(x2nbuff, 1, x_con_s_down_q);


	//aie::print(x_con_f, true, "_con_f ");

	//aie::vector<int32,16> x_vec = x2buff;

	aie::vector<int32,16> x_vec = x2buff - x2nbuff;
	//aie::print(x_vec, true, "tedt ");





	// 4-point NTT stage01

	v16int32 temp_buff = aie::transpose(x_vec, 4, 4) ;



	v4int32 k0, k1, k2, k3;

	k0 = srs(lmul4_sym     (temp_buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 0) , 0) ; // *yy++ = k0;
	k1 = srs(lmul4_antisym (temp_buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 1) , 0) + p2 ; // *yy++ = (k1+p2);

	k2 = srs(lmul4_sym     (temp_buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 2)  ,0) + p1 ; // *yy++ = (k2+p1);
	k3 = srs(lmul4_antisym (temp_buff, 0, 0x3210, 4, 12, zbuff, 0, 0x0000, 3)  ,0) + p2 ; // *yy   = (k3+p2);

	aie::vector<int32,8> x_vec2_first;
	x_vec2_first =   upd_v(x_vec2_first, 0, k0);
	x_vec2_first =   upd_v(x_vec2_first, 1, k1);

	aie::vector<int32,8> x_con_f2 = (srs(aie::mul(x_vec2_first, 255),0)); // mul saved in acc80, v16acc80 not possible; no srs() for v16int32
	v8int32 x_con_f2_down_q = srs( mul(aie::downshift(x_con_f2,16), p0), 0);

	aie::vector<int32,8> x_vec2 = x_vec2_first - x_con_f2_down_q;
	*yy++ = x_vec2;


	aie::vector<int32,8> x_vec2_second;

	x_vec2_second =  upd_v(x_vec2_second, 0, k2);
	x_vec2_second =  upd_v(x_vec2_second, 1, k3);

	aie::vector<int32,8> x_con_s2 = (srs(aie::mul(x_vec2_second, 255),0));
	v8int32 x_con_s2_down_q = srs( mul(aie::downshift(x_con_s2,16), p0), 0);

	aie::vector<int32,8> x_vec4 = x_vec2_second - x_con_s2_down_q;
	*yy = x_vec4;

}

#endif
