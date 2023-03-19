//#####################################################################
// Main
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#include <iostream>
#include "GridFluid.h"
#include "GridFluidDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc, char* argv[])
{
	GridFluidDriver driver;
	driver.Initialize();
	driver.Run();
}

#endif