#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "convection.h"
#include "domain.h"


void
CONVECTION_moveParticles( double **position, double **velocity){

  int iParticle;
  int iDim;

#pragma omp parallel for private(iDim)
  for( iParticle=0; iParticle < particle.totalNumber; iParticle++){
    if( particle.type[iParticle] == GHOST ) continue;

    /*
    if(parameter.flagOfTankCalculation == ON){
      if( particle.type[iParticle] == parameter.dummyWallType ){
	continue;
      }
    }
    */

    /*
    if(parameter.flagOfTankCalculation == ON){
      if( (particle.type[iParticle] == parameter.wallType )||(particle.type[iParticle] == parameter.dummyWallType )){
	if( (particle.flagOfFixVelocity[iParticle] == OFF )&&( particle.flagOfFixVelocityAndPressure_forNonReflectingBoundary[iParticle] == OFF)){
	    continue;
	}
      }
    }
    */

    for( iDim=0; iDim < NumberOfDimensions; iDim++){
      position[iDim][iParticle] +=  (velocity[iDim][iParticle] * timer.dt);
    }

  }


}


