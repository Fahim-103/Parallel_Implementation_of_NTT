#ifndef __NTTSUB_H__
#define __NTTSUB_H__


#include <adf.h>
//#include "twiddle.h"

#include "aie_api/aie.hpp"
#include "aie_api/utils.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/operators.hpp"


//#define rnd_floor 0


inline __attribute__((always_inline)) v8int32 modular_ntt ( v8int32 v){

	//modulus sub generation for mod = t - q [ (65535 * t) >> 32 ]

	int32 rr[8] = {65535, 65535, 65535, 65535, 65535, 65535, 65535, 65535};  //w'/q
	v8int32  r = *((v8int32*)rr);

	int32 po[8] = {65537, 65537, 65537, 65537, 65537, 65537, 65537, 65537};  // q
	v8int32 p0 = *((v8int32*)po);  // q

	v8int32 vcon, vr;
	vcon = srs(mul(srs(mul(v,r),32),p0),0) ;

	vr   = v - vcon;

	return vr;

}

inline __attribute__((always_inline)) void stage0_ntt ( int32 * restrict x,  int32 * restrict y,
														int32 * restrict tw1, int32 * restrict tw2, int32 * restrict tw3){

	unsigned int a =  get_srs_sat();
	printf("mode sat: %u \n", a);


	v8int32 * xx = (v8int32 *) x; //16 byte 4 elem update
	v8int32 * yy = (v8int32 *) y;

	v8int32 * t1 = (v8int32 *) tw1;
	v8int32 * t2 = (v8int32 *) tw2;
	v8int32 * t3 = (v8int32 *) tw3;

	v32int32 x1buff;

	// 4-point NTT twiddle
	int32 const chess_storage(%chess_alignof(v4int32)) tw[8] =   {1, 256, -1, -256,
											                      0,  0,   0,  0}; //for zbuff twiddle for 4-point NTT (w_4)^1 = (w_16)^4 = 16
	v8int32 * zz = (v8int32 *)tw; ////////
	v8int32 zbuff = *zz;

	//modular pre-requirements
	int32 pa[8] = {131074, 131074, 131074, 131074, 131074, 131074, 131074, 131074};  // 2q
	int32 pb[8] = {16843009, 16843009, 16843009, 16843009, 16843009, 16843009, 16843009, 16843009}; // q(w+1) ; w is from zbuff
	v8int32 p1 = *((v8int32*)pa);  // 2q
	v8int32 p2 = *((v8int32*)pb);  // q(w+1) ; w is from zbuff


	v8int32 m0, m1, m2, m3, m1_test, m2_test, m3_test, mr0, mr1, mr2, mr3, mrt, mrtf;

	v4int32 m00, m10, m20, m30, m01, m11, m21, m31;

	aie::vector<int32,8> vtest, v, vt, vrt, m1_s, m2_s, m3_s, mr1f_s, mr2f_s, mr3f_s;

	aie::mask<8> msk_lt;
	aie::vector<int32,8> c1=aie::broadcast<int32,8>(1);

	for (int i = 0; i < 32; i++){

		x1buff = upd_w(x1buff, 0, *xx); xx += 32;
		x1buff = upd_w(x1buff, 1, *xx); xx += 32;
		x1buff = upd_w(x1buff, 3, *xx); xx += 32;
		x1buff = upd_w(x1buff, 2, *xx); xx -= 95;
		//aie::vector<int32,32> x_vec = x1buff;
		//aie::print(x_vec, true, "x_vec ");

		// block of first second third fourt points, with tw multiplication


		m00 =           srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 0) , 0);
		m01 =           srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 0) , 0);


		m10 = srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 1) , 0) ; //first points of 4-blocks 4-point
		m11 = srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 1) , 0) ;


		m20 = srs(lmul4_sym     (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 2)  ,0) ;
		m21 = srs(lmul4_sym     (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 2)  ,0) ;


		m30 = srs(lmul4_antisym (x1buff, 0, 0x3210, 8, 24, zbuff, 0, 0x0000, 3)  ,0) ;
		m31 = srs(lmul4_antisym (x1buff, 4, 0x3210, 8, 28, zbuff, 0, 0x0000, 3)  ,0) ;


		//===========================================================================================================



		m0 = upd_v(m0, 0, m00);
		m0 = upd_v(m0, 1, m01);

		m1 = upd_v(m1, 0, m10);
		m1 = upd_v(m1, 1, m11);

		m2 = upd_v(m2, 0, m20);
		m2 = upd_v(m2, 1, m21);

		m3 = upd_v(m3, 0, m30);
		m3 = upd_v(m3, 1, m31);



		//===========================================================================================================

		m1_test = m1 + p2;
		m2_test = m2 + p1;
		m3_test = m3 + p2;

		//===========================================================================================================

		v = m1;
		vt = m1_test;
		msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
		m1_s=aie::select(v,vt,msk_lt); // select v2[i] if msk_lt[i] is true; select v1[i] if msk_lt[i] is	false;
		//aie::print(m1_s, true, "m1s " );

		//===========================================================================================================

		v = m2;
		vt = m2_test;
		msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
		m2_s=aie::select(v,vt,msk_lt); // select v2[i] if msk_lt[i] is true; select v1[i] if msk_lt[i] is	false;
		//aie::print(m2_s, true, "m1s " );

		//===========================================================================================================

		v = m3;
		vt = m3_test;
		msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
		m3_s=aie::select(v,vt,msk_lt); // select v2[i] if msk_lt[i] is true; select v1[i] if msk_lt[i] is	false;
		//aie::print(m3_s, true, "m1s " );

		mr0 = modular_ntt(m0); *yy = mr0; yy +=32;

		mr1 = srs(mul(modular_ntt(m1_s), *t1++),0)  ;
		mr2 = srs(mul(modular_ntt(m2_s), *t2++),0)  ;
		mr3 = srs(mul(modular_ntt(m3_s), *t3++),0)  ;

		//===========================================================================================================


		mrt  = (modular_ntt(mr1));
		mrtf = mrt+c1;

		v = mr1; ////////
		vt = mrt;
		vrt = mrtf;
		msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
		mr1f_s=aie::select(vt,vrt,msk_lt); // select v2[i] if msk_lt[i] is true; select v1[i] if msk_lt[i] is	false;
		//aie::print(mr1f_s, true, "mr1f_s " );
		*yy = mr1f_s; yy +=32;

		//===========================================================================================================


		mrt  = (modular_ntt(mr2)); ////////
		mrtf = mrt+c1;

		v = mr2; ////////
		vt = mrt;
		vrt = mrtf;
		msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
		mr2f_s=aie::select(vt,vrt,msk_lt); // select v2[i] if msk_lt[i] is true; select v1[i] if msk_lt[i] is	false;
		//aie::print(mr2f_s, true, "mr2f_s " );
		*yy = mr2f_s; yy +=32;

		//===========================================================================================================


		mrt  = (modular_ntt(mr3)); ////////
		mrtf = mrt+c1;

		v = mr3; ////////
		vt = mrt;
		vrt = mrtf;
		msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
		mr3f_s=aie::select(vt,vrt,msk_lt); // select v2[i] if msk_lt[i] is true; select v1[i] if msk_lt[i] is	false;
		//aie::print(mr3f_s, true, "mr2f_s " );
		*yy = mr3f_s; yy -=95;

}


}





#endif
