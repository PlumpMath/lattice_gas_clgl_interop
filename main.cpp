#include "Lattice_gas.h"


int main(){

	std::cout << "lattice gas simulation" << std::endl;

	//cpu
	std::cout << "cpu simulation" << std::endl;
	simul_cpu();

	//gpu
	std::cout << "gpu simulation" << std::endl;
	simul_gpu();	

	std::cin.get();

	return 0;
}