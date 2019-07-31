#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "file.h"
#include "copy.h"
#include "bucket.h"
#include "domain.h"
#include "neigh.h"
#include "density.h"
#include "timer.h"
#include "other.h"
#include "gravity.h"
#include "convection.h"
#include "collision.h"
#include "pressure.h"
#include "gradient.h"
#include "viscosity.h"
#include "sort.h"
#include "object.h"
#include "forcedMotion.h"
#include "memory.h"
#include "nzero.h"
#include "distance.h"
#include "maxmin.h"
#include "init.h"
#include "bubble.h"
#include "concentration.h"
#include "interface.h"
#include "vortex.h"


void
INIT_initializeParameters( int argumentCount,char **argumentVector ){

  timer.clockAtStartingTime = clock();
  timer.realTimeStart = TIMER_getWTime();

  timer.iTimeStep = 0;
  timer.iTimeStep_copy = -1;

  parameter.threads = 4;

  FILE_readInputFiles( argumentCount, argumentVector );

  INIT_initializeArguments( argumentCount, argumentVector );

  /*--- Memory of structParticle is allocated in the above FILE_readInputFiles(). ---*/

  COPY_storeInitialParticleProperty();

  DISTANCE_transformUnitOfDistanceFromRatioToMeter();

  MAXMIN_setMaxMinOfRadii();

  TIMER_initializeTimeFunctions();

  NZERO_calculateNZeroAndLambda( particle.averageDistance );

  if(parameter.flagOfBubbleCalculation==ON){
     Bubble_setBetaZero();
  }

  if (parameter.flagOfConcentrationCalculation == ON){
	  Concentration_initializeMixingIndexFile();
  }

  if (parameter.flagOfInterfaceCalculation == ON)
	  INTERFACE_initializeInterfaceFile();

  NEIGH_initializeNeighborTable();

  MEMORY_allocateMemoryForCoefficientMatrixOfPressure();

  if(parameter.flagOfForcedMotionOfRigidBody == ON ){
    FORCEDMOTION_initializeForcedMotion(  &forcedMotionOfRigidBody
					  ,parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody, parameter.typeNumberOfRigidParticle_forForcedMotion );

    /*
    FORCEDMOTION_initializeForcedMotion(  &forcedMotionOfRigidBody
					  ,parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody );
    */
    FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( &forcedMotionOfRigidBody );
  }

  COPY_setInitialCoodinatesAndVelocities();


  DOMAIN_initializeDomain();

  BUCKET_initializeBucket();

  NEIGH_setNeighborTable( particle.position );

  DENSITY_calculateParticleNumberDensity( particle.position );

  FILE_countNumberOfParticlesEachType();

  FILE_displayNumberOfParticles();

  fflush(FpForLog);

  //FILE_writeProfFile();

  NEIGH_setNeighborTable( particle.position );

  timer.clockAtBeginningOfMainLoop = clock();
  timer.realTimeLoop = TIMER_getWTime();

  timer.iTimeStep_copy++;

  if (parameter.flagOfConcentrationCalculation == ON){
	  Concentration_averageConcentration();
  }
  
  if (timer.simulationTime == 0)
	INTERFACE_CalculateInterfaceArea();

  VORTEX_setVelocityOfVortex();

  VORTEX_redistributeParticles();

  FILE_writeProfFile();

}

void
INIT_initializeArguments( int argumentCount,char **argumentVector ){

  int i;

  parameter.flagOfHideOutput = OFF;
  parameter.threads = 1;

  for(i = 1; i < argumentCount; i++ ) {
        if (strcmp( argumentVector[i], "-th" ) == 0 ) {
            parameter.threads = atoi(argumentVector[i+1]);
            i++;
        }
        if (strcmp( argumentVector[i], "-hide" ) == 0 ) {
            parameter.flagOfHideOutput = atoi(argumentVector[i+1]);
            i++;
        }
    }

}
