#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "struct.h"
#include "file.h"
#include "geometry.h"


# define PI           3.14159265358979323846  /* pi */

using namespace std;

void
INIT_initializeParameters(int argc,char **argv){

    strcpy(parameter.nameOfInputFile, "inputArguments.data");

	parameter.eps = 0.00001;

	parameter.width = 0.5;

	parameter.length = 2.0;

	parameter.max_num_particles = 2000000;

    FILE_readDataFile();

	FILE_showInput();

	GEOMETRY_initialize();
}

