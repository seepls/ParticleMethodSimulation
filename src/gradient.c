#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "pressure.h"
#include "memory.h"
#include "weight.h"
#include "mathe.h"
#include "convection.h"
#include "neigh.h"
#include "other.h"
#include "gradient.h"



void
GRADIENT_correctParticleVelocityAndPositionUsingPressureGradient( void ){

  OTHER_fill2dimDoubleArrayWithZero( NumberOfDimensions, particle.totalNumber, particle.velocity_correction);

  GRADIENT_calculatePressureGradientAndVelocityCorrection();


  MATHE_add2dimArrayTo2dimArray(  NumberOfDimensions, particle.totalNumber
				 ,particle.velocity, particle.velocity, particle.velocity_correction );

  CONVECTION_moveParticles( particle.position, particle.velocity_correction );

}



void
GRADIENT_calculatePressureGradientAndVelocityCorrection( void ){

  int iParticle, jParticle, iNeigh;

  double xji;
  double yji;
  double zji = 0.0;
  double d_vx;
  double d_vy;
  double d_vz = 0.0;

  double distanceIJ;
  double distanceIJ_squared;
  double absoluteValueOfVelocityCorrection;
  double averageDensity;

  int    *numberOfNeighborParticles;
  int    **neighborTable;
  int    **neighborTablePeriodic;


  NEIGH_selectNeighborTable(   &numberOfNeighborParticles
			       ,&neighborTable
			       ,parameter.radiusOfGradient_ratio
				   ,&neighborTablePeriodic
			       );

#pragma omp parallel for private(jParticle, iNeigh, xji, yji, zji, d_vx, d_vy, d_vz, distanceIJ, distanceIJ_squared, absoluteValueOfVelocityCorrection, averageDensity)
  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    particle.velocity_correction[XDIM][iParticle] = 0.0;
    particle.velocity_correction[YDIM][iParticle] = 0.0;

    if(NumberOfDimensions == 3){
      particle.velocity_correction[ZDIM][iParticle] = 0.0;
    }


    if( particle.flagOfBoundaryCondition[iParticle] == GHOST_OR_DUMMY )continue;
    if( particle.type[iParticle] == GHOST                   )continue;
    if( particle.type[iParticle] == parameter.dummyWallType )continue;
    if( particle.type[iParticle] == parameter.wallType      )continue;

    absoluteValueOfVelocityCorrection = 0.0;

    for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
      jParticle = neighborTable[iParticle][iNeigh];

      if(particle.type[jParticle]==GHOST || particle.type[jParticle]== parameter.dummyWallType )continue;

      xji = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
	  if(neighborTablePeriodic[iParticle][iNeigh]==1)
		  xji += parameter.pipeLength;
	  else if(neighborTablePeriodic[iParticle][iNeigh]==2)
		  xji -= parameter.pipeLength;
      yji = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);


      distanceIJ_squared = xji*xji + yji*yji;
      
      if(NumberOfDimensions == 3){
		zji = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
		distanceIJ_squared += zji*zji;
      }


      distanceIJ = sqrt(distanceIJ_squared);
      averageDensity = MATHE_average( physicalProperty.massDensity[particle.type[iParticle]]
				      ,physicalProperty.massDensity[particle.type[jParticle]] ); 

      absoluteValueOfVelocityCorrection  =  NumberOfDimensions * timer.dt;
      //absoluteValueOfVelocityCorrection *= (particle.pressure[jParticle] - particle.minPressureAmongNeighbors[iParticle])/distanceIJ;
      /*new*/
      absoluteValueOfVelocityCorrection *= (particle.pressure[jParticle] + particle.pressure[iParticle])/distanceIJ;


      absoluteValueOfVelocityCorrection *= WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfGradient);


      absoluteValueOfVelocityCorrection /= (averageDensity * parameter.nZeroOfGradient);


      d_vx = absoluteValueOfVelocityCorrection * xji / distanceIJ;
      d_vy = absoluteValueOfVelocityCorrection * yji / distanceIJ;

      if(NumberOfDimensions == 3){
	d_vz = absoluteValueOfVelocityCorrection * zji / distanceIJ;
      }

      particle.velocity_correction[XDIM][iParticle] += d_vx;
      particle.velocity_correction[YDIM][iParticle] += d_vy;

      if(NumberOfDimensions == 3){
	particle.velocity_correction[ZDIM][iParticle] += d_vz;
      }


    }


  }
}

