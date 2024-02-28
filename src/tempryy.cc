

v4int32 mlt, mlt_mod, v, vt;
aie::mask<4> msk_lt;
aie::vector<int32,4> c1=aie::broadcast<int32,8>(1);

mlt = srs(mul(m1, t1),0);
mlt_mod = modular_ntt_4(mlt);
//vmsk = mlt;
v  = mlt_mod;
vt = v + c1;
msk_lt=aie::lt(mlt,(int32)0);  //compare each element, true ifv1[i]<v2[i]
mout=aie::select(v,vt,msk_lt);

return mout;


inline __attribute__((always_inline)) v4int32 twdl_mod_mul_4 ( v4int32 m1, v4int32 t1){

	v4int32 mlt, mlt_mod, v, vt;
	aie::mask<4> msk_lt;
	aie::vector<int32,4> c1=aie::broadcast<int32,8>(1);

	mlt = srs(mul(m1, t1),0);
	mlt_mod = modular_ntt_4(mlt);
	//vmsk = mlt;
	v  = mlt_mod;
	vt = v + c1;
	msk_lt=aie::lt(mlt,(int32)0);  //compare each element, true ifv1[i]<v2[i]
	mout=aie::select(v,vt,msk_lt);

	return mout;
}

inline __attribute__((always_inline)) v4int32 twdl_mod_sub_4 ( v4int32 m, v4int32 p){

	v4int32  v, vt;
	aie::mask<4> msk_lt;

	v = modular_ntt_4(m);
	vt = v + p;
	msk_lt=aie::lt(v,(int32)0);  //compare each element, true ifv1[i]<v2[i]
	mout=aie::select(v,vt,msk_lt);

	return mout;
}
