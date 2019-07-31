#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "gravity.h"



void
	GRAVITY_calculateGravity( void ){

		int iParticle;
		int iDim;

#pragma omp parallel for private(iDim)
		for(iParticle= 0; iParticle < particle.totalNumber; iParticle++){

			if ( particle.type[iParticle]==GHOST ) continue;
			if ( particle.type[iParticle]==parameter.dummyWallType ) continue;
			if ( particle.type[iParticle]==parameter.wallType ) continue;

			for(iDim= 0; iDim < NumberOfDimensions; iDim++){
				particle.velocity[iDim][iParticle] += physicalProperty.gravity[iDim] * timer.dt;
			}
		}

}
