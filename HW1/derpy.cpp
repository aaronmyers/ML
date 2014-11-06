#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>

int main(int argc,char* argv[]) {
	if(argc!=5) {
		std::cout << "Ya messed up, boy" << std::endl;

	}
	int rank=argv[1];
	double lambda=argv[2];
	int nr_threads=argv[3];
	std::string data_dir=argv[4];

	




	return 0;
}
