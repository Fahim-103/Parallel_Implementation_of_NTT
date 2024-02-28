#include "graph.h"

class TEST: public graph {

public:

	input_plio  plin[4];
	output_plio plout[4];

	//graph object instant
	ntt4k_graph g;



    TEST(){

		plin[0] = input_plio::create("in00a",  plio_32_bits, "data/in00a.txt", 500);
		plin[1] = input_plio::create("in00b",  plio_32_bits, "data/in00b.txt", 500);
		plin[2] = input_plio::create("in00c",  plio_32_bits, "data/in00c.txt", 500);
		plin[3] = input_plio::create("in00d",  plio_32_bits, "data/in00d.txt", 500);

		plout[0] = output_plio::create("out00a",  plio_32_bits, "data/out00a.txt", 500);
		plout[1] = output_plio::create("out00b",  plio_32_bits, "data/out00b.txt", 500);
		plout[2] = output_plio::create("out00c",  plio_32_bits, "data/out00c.txt", 500);
		plout[3] = output_plio::create("out00d",  plio_32_bits, "data/out00d.txt", 500);




		connect<>(plin[0].out[0], g.in[0]);
		connect<>(plin[1].out[0], g.in[1]);
		connect<>(plin[2].out[0], g.in[2]);
		connect<>(plin[3].out[0], g.in[3]);


		connect<>(g.out[0], plout[0].in[0]);
		connect<>(g.out[1], plout[1].in[0]);
		connect<>(g.out[2], plout[2].in[0]);
		connect<>(g.out[3], plout[3].in[0]);



		}

};
