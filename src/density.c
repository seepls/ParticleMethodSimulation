#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "domain.h"
#include "weight.h"
#include "distance.h"
#include "maxmin.h"
#include "neigh.h"
#include "density.h"


void
DENSITY_calculateParticleNumberDensity( double **position ){

  int iParticle, jParticle;

  int iType;
  int iNeigh;

  double distanceIJ;
  double weightIJ;

  int *numberOfNeighborParticles;
  int **neighborTable;
  int **neighborTablePeriodic;

  /*
  numberOfNeighborParticles = NULL;
  neighborTable             = NULL;
  */
  /*
  numberOfNeighborParticles = particle.numberOfNeighborParticles_large;
  neighborTable             = particle.neighborTable_large;
  */


  NEIGH_selectNeighborTable( &numberOfNeighborParticles
							 ,&neighborTable
							 ,parameter.radiusOfParticleNumberDensity_ratio
							 ,&neighborTablePeriodic
							 );

  /*
  NEIGH_selectNeighborTable(  numberOfNeighborParticles
							  ,neighborTable
							  ,parameter.radiusOfParticleNumberDensity_ratio
							  );
  */

  for(iParticle=0; iParticle <particle.totalNumber; iParticle++){

    if(particle.type[iParticle] == GHOST) continue;

    particle.particleNumberDensity[iParticle] = 0.0;
    particle.particleNumberDensity_withConstantWeight[iParticle] = 0.0;

    for(iType=0; iType < parameter.numberOfParticleTypes; iType++) {
      particle.particleNumberDensity_eachParticleType[iParticle][iType] = 0.0;
    }

    for(iNeigh=0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
      jParticle  = neighborTable[iParticle][iNeigh];
      //distanceIJ = DISTANCE_calculateDistanceBetweenParticles( position, iParticle, jParticle );
      distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic( position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);
      weightIJ   = WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfParticleNumberDensity );
      particle.particleNumberDensity_eachParticleType[iParticle][particle.type[jParticle]] += weightIJ;

	  /*if((iParticle==0)||(iParticle==2000))
		  printf("iParticle: %d, jParticle: %d, iNeigh: %d BoundaryNeighborTable: %d, distanceIJ: %f, weightIJ: %f particle.particleNumberDensity_eachType: %f \n", 
		  iParticle,jParticle,iNeigh,neighborTablePeriodic[iParticle][iNeigh],distanceIJ,weightIJ,particle.particleNumberDensity_eachParticleType);*/

      if( distanceIJ <= 2.1 * particle.averageDistance){
	 particle.particleNumberDensity_withConstantWeight[iParticle] += 1.0;
      }
    }

    for(iType=0; iType < parameter.numberOfParticleTypes; iType++) {
      particle.particleNumberDensity[iParticle] += particle.particleNumberDensity_eachParticleType[iParticle][iType];
    }

  }


}


/*
void
DENSITY_setStandardParticleNumberDensity( void ){

  MAXMIN_updateMaxMinOfParticleProperties();

  parameter.nZeroOfParticleNumberDensity = particle.maxParticleNumberDensity;
  parameter.nZeroOfGradient              = particle.maxParticleNumberDensity;

  DENSITY_displayStandardNumberDensity();

}
*/


/*
void
DENSITY_calculateLambdaOfLaplacianOperator( void ){

  int    iParticle,jParticle;
  int    iNeigh;
  double distanceIJ, distanceIJ_squared;
  double weightIJ;
  double totalWeight;
  double tatalWeightTimesSquaredDistance;

  iParticle =particle.standardParticleNumber;

  totalWeight    = 0.0;
  tatalWeightTimesSquaredDistance = 0.0;

  int    *numberOfNeighborParticles;
  int    **neighborTable;

  numberOfNeighborParticles = particle.numberOfNeighborParticles_large;
  neighborTable           = particle.neighborTable_large;

  for(iNeigh=0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
    jParticle  = neighborTable[iParticle][iNeigh];

    if(particle.type[jParticle] != GHOST ){

      distanceIJ_squared = DISTANCE_calculateSquaredDistance(particle.position, iParticle, jParticle );
      distanceIJ  = sqrt(distanceIJ_squared);

      weightIJ   = WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure );
      tatalWeightTimesSquaredDistance += weightIJ * distanceIJ_squared;
      totalWeight    += weightIJ;

    }

  }

  parameter.nZeroOfLaplacianForPressure = totalWeight;

  parameter.lambda             = tatalWeightTimesSquaredDistance / totalWeight;
  parameter.lambdaOfLaplacianForPressure  = parameter.lambda * parameter.relaxationCoefficientOfLambda;
  parameter.lambdaOfLaplacianForViscosity = parameter.lambda;
  parameter.lambdaTimesNZeroForViscosity  = tatalWeightTimesSquaredDistance;

  DENSITY_displayCalculatedLambda();

}
*/




void
DENSITY_displayCalculatedLambda( void ){

  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        Lambda of Laplacian operator                \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"lambdaTimesNZeroForViscosity  = %.10e\n", parameter.lambdaTimesNZeroForViscosity );
  fprintf(FpForLog,"lambdaOfLaplacianForPressure  = %.10e\n", parameter.lambdaOfLaplacianForPressure  );
  fprintf(FpForLog,"lambdaOfLaplacianForViscosity = %.10e\n", parameter.lambdaOfLaplacianForViscosity );
  fprintf(FpForLog,"relaxationCoefficientOfLambda   = %.10e\n", parameter.relaxationCoefficientOfLambda   );
  fprintf(FpForLog,"standard particle             = %d th particle\n",particle.standardParticleNumber );
  fprintf(FpForLog,"nZeroOfLaplacianForPressure   = %.10e\n", parameter.nZeroOfLaplacianForPressure  );
         
  fprintf(FpForLog,"\n\n");

}


void
DENSITY_displayStandardNumberDensity( void ){

  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        n0: standard number density                 \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"parameter.nZeroOfParticleNumberDensity =  %lf [particles]\n", parameter.nZeroOfParticleNumberDensity );
  fprintf(FpForLog,"parameter.nZeroOfGradient              =  %lf [particles]\n", parameter.nZeroOfGradient );

  fprintf(FpForLog,"\n\n");

}


