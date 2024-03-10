#include <test.h>

TEST dut;

int main(void) {
	dut.init();
    dut.run(1) ;
	dut.end() ;
    return 0 ;
}
