#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "neigh.h"
#include "density.h"
#include "distance.h"
#include "weight.h"
#include "mathe.h"
#include "solver.h"
#include "domain.h"
#include "copy.h"
#include "neigh.h"
#include "memory.h"
#include "other.h"
#include "pressure.h"



void
PRESSURE_calculatePressure( void ){

  int flagOfPressureCalculation;

  NEIGH_setNeighborTable( particle.position );

  DENSITY_calculateParticleNumberDensity( particle.position );

  PRESSURE_setDirichletBoundaryCondition( );


  flagOfPressureCalculation = PRESSURE_checkNecessityOfPressureCalculation( );

  if( flagOfPressureCalculation == ON ){

    PRESSURE_setSourceTerm();

    PRESSURE_setCoefficientMatrixOfPressure();

    SOLVER_solveSimultaniousEquations();

    if(parameter.flagOfHighViscosityCalculation == ON){
       PRESSURE_correctPressure();
    }

    PRESSURE_setMinusPressureZero();

    PRESSURE_setMinimumPressureAmongNeighbors();

  }else{

    OTHER_fillOneDimDoubleArrayWithZero( particle.totalNumber, particle.pressure );

    fprintf(FpForLog,"WARNING: Pressure was not calculated because all particles were surface-particles or wall-particles.\n");

    PRESSURE_setMinimumPressureAmongNeighbors();

  }

}




void
PRESSURE_setSourceTerm( void ){

  int iParticle;

  double n0         = parameter.nZeroOfParticleNumberDensity;
  double dt_squared = timer.dt * timer.dt;

  /*after*/
  int    jParticle;
  int    iNeigh;
  int    iDim;

  double distanceIJ;
  double distanceIJ_squared;
  double sigma;

  double dt = timer.dt;
  double weightIJ;

  int    *numberOfNeighborParticles;
  int    **neighborTable;
  int    **neighborTablePeriodic;

  NEIGH_selectNeighborTable(  &numberOfNeighborParticles
                              ,&neighborTable
                              ,parameter.radiusOfLaplacianForPressure_ratio
							  ,&neighborTablePeriodic
                              );


  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if ( (particle.type[iParticle]==GHOST) || (particle.type[iParticle] == parameter.dummyWallType) ){

      particle.sourceTermOfPressure[iParticle] = 0.0;


    }else if( particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE){

        if(parameter.flagOfTanakaAndMasunagaModel==ON){

           sigma=0.0;

           for(iNeigh=0; iNeigh< numberOfNeighborParticles[iParticle];iNeigh++){

               jParticle=neighborTable[iParticle][iNeigh];

               if(particle.type[jParticle]==GHOST) continue;

               //distanceIJ = DISTANCE_calculateDistanceBetweenParticles( particle.position, iParticle, jParticle  );
			   distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);

               //distanceIJ_squared = DISTANCE_calculateSquaredDistance(particle.position, iParticle, jParticle);
			   distanceIJ_squared = DISTANCE_calculateSquaredDistancePeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);

               weightIJ = WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure);

               for(iDim= 0; iDim < NumberOfDimensions; iDim++){
                   sigma -= ((particle.velocity[iDim][jParticle]-particle.velocity[iDim][iParticle])*(particle.position[iDim][jParticle]-particle.position[iDim][iParticle])*weightIJ/distanceIJ_squared);
               }
           }
           particle.sourceTermOfPressure[iParticle]=(NumberOfDimensions*sigma/(dt*n0))+parameter.valueOfGamma*((particle.particleNumberDensity_previous[iParticle]-n0)/(n0*dt_squared));

        }else{
           particle.sourceTermOfPressure[iParticle] = ( particle.particleNumberDensity[iParticle] - n0 )/( dt_squared * n0 );
        }

    }else if(particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE){

      particle.sourceTermOfPressure[iParticle] = 0.0;

    }
  }
}





void
PRESSURE_setCoefficientMatrixOfPressure( void ){

  int    iParticle, jParticle;
  int    iNeigh;

  double n0     = parameter.nZeroOfLaplacianForPressure;
  double lambda = parameter.lambdaOfLaplacianForPressure;

  double distanceIJ;
  double coefficientOfParticleJ;
  double averageMassDensity;

  double dt_squared = (timer.dt * timer.dt);
  double weightIJ;


  int    *numberOfNeighborParticles;
  int    **neighborTable;
  int    **neighborTablePeriodic;


  NEIGH_selectNeighborTable(  &numberOfNeighborParticles
			      ,&neighborTable
			      ,parameter.radiusOfLaplacianForPressure_ratio
				  ,&neighborTablePeriodic
			      );

  for( iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.flagOfBoundaryCondition[iParticle] != INNER_PARTICLE ) continue;

    particle.coefficientMatrixOfPressure[iParticle][0]=0.0;

    for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
      jParticle = neighborTable[iParticle][iNeigh];

      if( particle.flagOfBoundaryCondition[jParticle] == GHOST_OR_DUMMY ){
	      particle.coefficientMatrixOfPressure[iParticle][iNeigh+1]=0.0;
      }else {
	     //distanceIJ = DISTANCE_calculateDistanceBetweenParticles( particle.position, iParticle, jParticle  );
		 distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);

	     averageMassDensity = MATHE_average( physicalProperty.massDensity[particle.type[iParticle]]
					    ,physicalProperty.massDensity[particle.type[jParticle]] );

	     weightIJ = WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure);

	     coefficientOfParticleJ  = 2.0 * NumberOfDimensions * weightIJ /( n0 * lambda  * averageMassDensity );

	     particle.coefficientMatrixOfPressure[iParticle][iNeigh+1]  = - coefficientOfParticleJ;

         if(parameter.flagOfTanakaAndMasunagaModel==ON){
            particle.coefficientMatrixOfPressure[iParticle][0]        +=   parameter.valueOfC*coefficientOfParticleJ;
         }else{
            particle.coefficientMatrixOfPressure[iParticle][0]        +=   coefficientOfParticleJ;
         }
      }
    }

    particle.coefficientMatrixOfPressure[iParticle][0] += physicalProperty.compressibility[particle.type[iParticle]]/( dt_squared );

  }

}



void
PRESSURE_setMinusPressureZero( void ){

  int iParticle;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.pressure[iParticle] < 0.0 ){
      particle.pressure[iParticle] = 0.0;
    }

  }

}



void
PRESSURE_setMinimumPressureAmongNeighbors( void ){

  int    iParticle,jParticle;
  int    iNeigh;

  double minimumPressure;
  double distanceIJ_squared;


  int    *numberOfNeighborParticles;
  int    **neighborTable;
  int    **neighborTablePeriodic;


  NEIGH_selectNeighborTable(   &numberOfNeighborParticles
							  ,&neighborTable
							  ,parameter.radiusOfGradient_ratio
							  ,&neighborTablePeriodic
							  );


  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.flagOfBoundaryCondition[iParticle] == GHOST_OR_DUMMY ) continue;
    if( particle.type[iParticle]        == parameter.wallType ) continue;

    minimumPressure = particle.pressure[iParticle];

    for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
      jParticle = neighborTable[iParticle][iNeigh];

      if(particle.type[jParticle] == GHOST              )continue;
      if(particle.type[jParticle] == parameter.dummyWallType)continue;

	  //distanceIJ_squared = DISTANCE_calculateSquaredDistance(particle.position, iParticle, jParticle);
	  distanceIJ_squared = DISTANCE_calculateSquaredDistancePeriodic(particle.position, iParticle, jParticle,neighborTablePeriodic[iParticle][iNeigh]);
	  if(distanceIJ_squared > parameter.radiusOfGradient_squared) continue;

      if(particle.pressure[jParticle] < minimumPressure) minimumPressure = particle.pressure[jParticle];

    }

    particle.minPressureAmongNeighbors[iParticle] = minimumPressure;

  }
}



int
PRESSURE_checkNecessityOfPressureCalculation( void ){

  int iParticle;
  int flagOfPressureCalculation = OFF;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE){
      flagOfPressureCalculation = ON;
    }

  }

  return flagOfPressureCalculation;

}





void
PRESSURE_setDirichletBoundaryCondition( void ){

  if(parameter.flagOfTanakaAndMasunagaModel==ON){
     PRESSURE_findParticlesBeingOnFreeSurfaceTheNumberOfNeighParticleBased();
  }else{
     PRESSURE_findParticlesBeingOnFreeSurface();
  }

  PRESSURE_checkThatDirichletBoundaryConditionIsConnected();

}





void
PRESSURE_findParticlesBeingOnFreeSurface( void ){

  int    iParticle;
  double n0 = parameter.nZeroOfParticleNumberDensity;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( (particle.type[iParticle] == GHOST) || (particle.type[iParticle] == parameter.dummyWallType)){

      particle.flagOfBoundaryCondition[iParticle] = GHOST_OR_DUMMY;

    }else if( (particle.particleNumberDensity[iParticle] / n0 ) < parameter.thresholdOfParticleNumberDensity_ratio ){

      particle.flagOfBoundaryCondition[iParticle] = SURFACE_PARTICLE;

    }else{

      particle.flagOfBoundaryCondition[iParticle] = INNER_PARTICLE;

    }
  }

}

void
PRESSURE_findParticlesBeingOnFreeSurfaceTheNumberOfNeighParticleBased( void ){

    int    iParticle;
    double N0 = parameter.nZeroOfnumberOfNeighborParticles;
    double numberOfNeighborParticles_ratio=0.80;

    int    *numberOfNeighborParticles;
    int    **neighborTable;
	int    **neighborTablePeriodic;

    NEIGH_selectNeighborTable(   &numberOfNeighborParticles
							  ,&neighborTable
							  ,parameter.radiusOfParticleNumberDensity_ratio
							  ,&neighborTablePeriodic
							  );


    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

        if( (particle.type[iParticle] == GHOST) || (particle.type[iParticle] == parameter.dummyWallType)){

            particle.flagOfBoundaryCondition[iParticle] = GHOST_OR_DUMMY;

        }else if( numberOfNeighborParticles[iParticle] < N0*numberOfNeighborParticles_ratio ){

            particle.flagOfBoundaryCondition[iParticle] = SURFACE_PARTICLE;

        }else{

            particle.flagOfBoundaryCondition[iParticle] = INNER_PARTICLE;

        }
    }

}



void
PRESSURE_checkThatDirichletBoundaryConditionIsConnected( void ){

  int iParticle, jParticle, iNeigh;
  int count;
  int *checkArray;
  char buf[256];
  int iLoop=0;

  int    *numberOfNeighborParticles = particle.numberOfNeighborParticles_large;
  int    **neighborTable            = particle.neighborTable_large;

  checkArray = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber,
	                   "checkArray [in function of checkThatDirichletBoundariesAreConnected()]");

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if (particle.flagOfBoundaryCondition[iParticle]== GHOST_OR_DUMMY ){

      checkArray[iParticle]= GHOST_OR_DUMMY;

    }else if (particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE){

      checkArray[iParticle] = CONNECTED;

    }else if (particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE){

      checkArray[iParticle] = NOT_CONNECTED;

    }else{

      sprintf(buf,"[in PRESSURE_checkThatDirichletBoundariesAreConnected]\n particle.surface is not adequate.\n particle.flagOfBoundaryCondition[%d] = %d", iParticle, particle.flagOfBoundaryCondition[iParticle]);

      OTHER_endProgram("buf");

    }

  }


  do {
    count=0;

    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

      if(checkArray[iParticle] == CONNECTED){

	for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
	  jParticle = neighborTable[iParticle][iNeigh];

	  if(checkArray[jParticle]== NOT_CONNECTED){
	    checkArray[jParticle] = CONNECTED;
	  }
	}

	checkArray[iParticle] = CHECKED;
	count ++;

      }
    }

    iLoop++;

  } while (count > 0);


  PRESSURE_ifBoundaryConditionIsNotConnected_increaseDiagonalElementsOfCoefficientMatrixOfPressureToSolveSimultaneousEquationsWithoutBoundaryCondition( checkArray );

  free(checkArray);

}


void
PRESSURE_ifBoundaryConditionIsNotConnected_increaseDiagonalElementsOfCoefficientMatrixOfPressureToSolveSimultaneousEquationsWithoutBoundaryCondition( int *checkArray ){

  int iParticle;

  double increaseRate = parameter.increaseRateOfDiagonalTermInCoefficientMatrixOfPressure;


  for( iParticle=0; iParticle < particle.totalNumber; iParticle++){
    if(checkArray[iParticle] == NOT_CONNECTED ){

      particle.coefficientMatrixOfPressure[iParticle][0] = increaseRate * particle.coefficientMatrixOfPressure[iParticle][0];

     // PRESSURE_displayWarnigMessageForNoDirichletBoundaryCondition( iParticle );
    }

  }


}



void
PRESSURE_displayWarnigMessageForNoDirichletBoundaryCondition( int iParticle ){

  fprintf(FpForLog,"WARNIG: No Dirichlet boundary condition.");
  fprintf(FpForLog,"  type[%d] = %d",iParticle, particle.type[iParticle]);

  fprintf(FpForLog,"  position[XDIM][%d] = %lf",iParticle, particle.position[XDIM][iParticle]);
  fprintf(FpForLog,"  position[YDIM][%d] = %lf",iParticle, particle.position[YDIM][iParticle]);

  if(NumberOfDimensions == 3){
    fprintf(FpForLog,"  position[ZDIM][%d] = %lf",iParticle, particle.position[ZDIM][iParticle]);
  }

  fprintf(FpForLog,"\n");

}

void
PRESSURE_correctPressure( void ){

    int iParticle;
    double kinematicViscosity = physicalProperty.kinematicViscosity;

    double n0         = parameter.nZeroOfParticleNumberDensity;
    double dt_squared = timer.dt * timer.dt;

    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

            particle.pressure[iParticle] += (1.0/2.0)*kinematicViscosity*timer.dt*( particle.particleNumberDensity[iParticle] - n0 )/( dt_squared * n0 );

    }
}
