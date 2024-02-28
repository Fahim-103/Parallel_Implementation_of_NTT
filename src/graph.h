#include <kernels.h>

extern int32 ntt_1a[1024];/////////////////////
extern int32 ntt_2a[1024];

extern int32 stage0_tw[768];
extern int32 stage1_tw[192];
extern int32 stage2_tw[48];



extern int32 last_stage_tw1[1024];
extern int32 last_stage_tw2[1024];
extern int32 last_stage_tw3[1024];


using namespace adf;


class sub1k_graph_a : public adf::graph{

public:

    parameter buf1;  ////////////////////////////////////
    parameter buf2;

    parameter buf3;
    parameter buf4;
    parameter buf5;

    //parameter buf6;

	adf::port<input>  in;
    adf::port<output> out;
    adf::kernel kl;

    sub1k_graph_a(){

    	kl = adf::kernel::create(ntt_a);
    	adf::runtime<ratio>(kl) = 0.8 ;
    	source(kl) = "ntt_a.cc" ;


        buf1 = parameter::array(ntt_1a);//////////////////////
        buf2 = parameter::array(ntt_2a);

        buf3 = parameter::array(stage0_tw);
        buf4 = parameter::array(stage1_tw);
        buf5 = parameter::array(stage2_tw);

       // buf6 = parameter::array(stage1_tw);


        connect<>(buf1, kl);///////////////////////////
        connect<>(buf2, kl);

        connect<>(buf3, kl);
        connect<>(buf4, kl);
        connect<>(buf5, kl);

        //connect<>(buf6, kl);



    	adf::connect<window<1024*4>>(in, kl.in[0]);
    	adf::connect<window<1024*4>>(kl.out[0], out);

    };

};


class sub1k_graph_b : public adf::graph{

public:

    parameter buf1;  ////////////////////////////////////
    parameter buf2;

    parameter buf3;
    parameter buf4;
    parameter buf5;

    //parameter buf6;

	adf::port<input>  in;
    adf::port<output> out;
    adf::kernel kl;

    sub1k_graph_b(){

    	kl = adf::kernel::create(ntt_a);
    	adf::runtime<ratio>(kl) = 0.8 ;
    	source(kl) = "ntt_b.cc" ;


        buf1 = parameter::array(ntt_1a);//////////////////////
        buf2 = parameter::array(ntt_2a);

        buf3 = parameter::array(stage0_tw);
        buf4 = parameter::array(stage1_tw);
        buf5 = parameter::array(stage2_tw);

       // buf6 = parameter::array(stage1_tw);


        connect<>(buf1, kl);///////////////////////////
        connect<>(buf2, kl);

        connect<>(buf3, kl);
        connect<>(buf4, kl);
        connect<>(buf5, kl);

        //connect<>(buf6, kl);



    	adf::connect<window<1024*4>>(in, kl.in[0]);
    	adf::connect<window<1024*4>>(kl.out[0], out);

    };

};


class sub1k_graph_c : public adf::graph{

public:

    parameter buf1;  ////////////////////////////////////
    parameter buf2;

    parameter buf3;
    parameter buf4;
    parameter buf5;

    //parameter buf6;

	adf::port<input>  in;
    adf::port<output> out;
    adf::kernel kl;

    sub1k_graph_c(){

    	kl = adf::kernel::create(ntt_a);
    	adf::runtime<ratio>(kl) = 0.8 ;
    	source(kl) = "ntt_c.cc" ;


        buf1 = parameter::array(ntt_1a);//////////////////////
        buf2 = parameter::array(ntt_2a);

        buf3 = parameter::array(stage0_tw);
        buf4 = parameter::array(stage1_tw);
        buf5 = parameter::array(stage1_tw);

        //buf6 = parameter::array(stage1_tw);


        connect<>(buf1, kl);///////////////////////////
        connect<>(buf2, kl);

        connect<>(buf3, kl);
        connect<>(buf4, kl);
        connect<>(buf5, kl);

        //connect<>(buf6, kl);



    	adf::connect<window<1024*4>>(in, kl.in[0]);
    	adf::connect<window<1024*4>>(kl.out[0], out);

    };

};

class sub1k_graph_d : public adf::graph{

public:

    parameter buf1;  ////////////////////////////////////
    parameter buf2;

    parameter buf3;
    parameter buf4;
    parameter buf5;

    //parameter buf6;

	adf::port<input>  in;
    adf::port<output> out;
    adf::kernel kl;

    sub1k_graph_d(){

    	kl = adf::kernel::create(ntt_a);
    	adf::runtime<ratio>(kl) = 0.8 ;
    	source(kl) = "ntt_d.cc" ;


        buf1 = parameter::array(ntt_1a);//////////////////////
        buf2 = parameter::array(ntt_2a);

        buf3 = parameter::array(stage0_tw);
        buf4 = parameter::array(stage1_tw);
        buf5 = parameter::array(stage2_tw);

        //buf6 = parameter::array(stage1_tw);


        connect<>(buf1, kl);///////////////////////////
        connect<>(buf2, kl);

        connect<>(buf3, kl);
        connect<>(buf4, kl);
        connect<>(buf5, kl);

        //connect<>(buf6, kl);



    	adf::connect<window<1024*4>>(in, kl.in[0]);
    	adf::connect<window<1024*4>>(kl.out[0], out);

    };

};

class last_graph : public adf::graph{

public:

    parameter buf1;  ////////////////////////////////////
    parameter buf2;
    parameter buf3;


	adf::port<input>  in[4];
    adf::port<output> out[4];
    adf::kernel kl;


    last_graph(){

    	kl = adf::kernel::create(ntt_last);
    	adf::runtime<ratio>(kl) = 0.8 ;
    	source(kl) = "ntt_last.cc" ;


        buf1 = parameter::array(last_stage_tw1);//////////////////////
        buf2 = parameter::array(last_stage_tw2);
        buf3 = parameter::array(last_stage_tw3);



        connect<>(buf1, kl);///////////////////////////
        connect<>(buf2, kl);
        connect<>(buf3, kl);



    	adf::connect<window<1024*4>>(in[0], kl.in[0]);
    	adf::connect<window<1024*4>>(in[1], kl.in[1]);
    	adf::connect<window<1024*4>>(in[2], kl.in[2]);
    	adf::connect<window<1024*4>>(in[3], kl.in[3]);


    	adf::connect<window<1024*4>>(kl.out[0], out[0]);
    	adf::connect<window<1024*4>>(kl.out[1], out[1]);
    	adf::connect<window<1024*4>>(kl.out[2], out[2]);
    	adf::connect<window<1024*4>>(kl.out[3], out[3]);

    };

};


class ntt4k_graph : public adf::graph{

public:



	adf::port<input>  in[4];
    adf::port<output> out[4];
    adf::kernel kl;

    sub1k_graph_a ka;	//graph_class_name   object_name
    sub1k_graph_b kb;
    sub1k_graph_c kc;
    sub1k_graph_d kd;
    last_graph    kz;

    ntt4k_graph(){



    	adf::connect<>(in[0], ka.in);
    	adf::connect<>(in[1], kb.in);
    	adf::connect<>(in[2], kc.in);
    	adf::connect<>(in[3], kd.in);

    	adf::connect<>(ka.out, kz.in[0]);
    	adf::connect<>(kb.out, kz.in[1]);
    	adf::connect<>(kc.out, kz.in[2]);
    	adf::connect<>(kd.out, kz.in[3]);

    	adf::connect<>(kz.out[0], out[0]);
    	adf::connect<>(kz.out[1], out[1]);
    	adf::connect<>(kz.out[2], out[2]);
    	adf::connect<>(kz.out[3], out[3]);


    };

};




