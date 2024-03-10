#include <adf.h>



using namespace adf;

void ntt_a (
		input_window_int32  * restrict inputx,
	    output_window_int32 * restrict outputy);

void ntt_b (
		input_window_int32  * restrict inputx,
	    output_window_int32 * restrict outputy);

void ntt_c (
		input_window_int32  * restrict inputx,
	    output_window_int32 * restrict outputy);

void ntt_d (
		input_window_int32  * restrict inputx,
	    output_window_int32 * restrict outputy);


void ntt_last (input_window_int32  * restrict inputx0,
			   input_window_int32  * restrict inputx1,
			   input_window_int32  * restrict inputx2,
			   input_window_int32  * restrict inputx3,

			   output_window_int32 * restrict outputy0,
			   output_window_int32 * restrict outputy1,
			   output_window_int32 * restrict outputy2,
			   output_window_int32 * restrict outputy3);

