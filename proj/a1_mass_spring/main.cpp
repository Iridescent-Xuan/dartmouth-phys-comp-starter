//#####################################################################
// Main
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#include <iostream>
#include <cstdio>
#include "MassSpringInteractiveDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc, char* argv[])
{
	////default arguments
	int test = 1;
	int scale = 1;

	////parse from command lines
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-test") == 0) {
			test = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-scale") == 0) {
			scale = atoi(argv[++i]);
		}
	}
	std::cout << "[Mass spring simulation driver arguments]: -test = " << test << ", -scale = " << scale << std::endl;

	////initialize driver
	MassSpringInteractivDriver driver;
	driver.scale = scale;
	driver.test = test;
	driver.Initialize();
	driver.Run();
}

#endif