#ifndef __NTT_LAST_SUB_H__
#define __NTT_LAST_SUB_H__

#include <adf.h>

#include "aie_api/aie.hpp"
#include "aie_api/utils.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/operators.hpp"


inline __attribute__((always_inline)) v8int32 twdl_mod_mul_8 ( v8int32 m, v8int32 t){   // lane by lane twidl multiplication

	aie::mask<8> msk_lt;
	aie::vector<int32,8> mlt, mlt_con, mlt_mod, mlt_mod_c1, r, p;;

	r=aie::broadcast<int32,8>(349496);
	p=aie::broadcast<int32,8>(12289);


	mlt = srs(mul(m, t),0);


	//modular_ntt //modulus sub generation for mod = m - q [ ( m * r ) >> 32 ]
	mlt_con = srs(mul(srs(mul(mlt,r),32),p),0) ;
	mlt_mod   = mlt - mlt_con;

	return mlt_mod;
}

inline __attribute__((always_inline)) v4int32 twdl_mod_sub_4 ( v4int32 m, v4int32 pn){   // checks if the radx-4 output is negative and performs mod

	aie::mask<4> msk_lt;
	aie::vector<int32,4> m_mod, m_p, msel, m_con, r, p, m_in;

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



#endif

