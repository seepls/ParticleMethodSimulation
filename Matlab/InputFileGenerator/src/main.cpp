//===================================================
//Program to create input files for Poisseuille Flow
//
//===================================================


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "struct.h"
#include "file.h"
#include "init.h"

#define PI 3.14159265

using namespace std;

int main(int argc, char** argv)
{
	//------------
	// get input
	//------------
	INIT_initializeParameters(argc,argv);

	//-------------------------
	// create the input files
	//-------------------------
	FILE_output_data_file();
	FILE_output_grid_bubble_file();
	if (!strcmp(parameter.flagOfForcedMotion,"on"))
		FILE_output_forcedMotion_file();

	cout << "Press Enter to Finish";
	cin.ignore();
}
