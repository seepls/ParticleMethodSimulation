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
#include "density.h"
#include "memory.h"
#include "other.h"
#include "nzero.h"


void
NZERO_calculateNZeroAndLambda( double averageDistanceBetweenParticles ){

  int     iParticle;
  int     jParticle;

  double  distanceIJ;
  double  distanceIJ_squared;

  double  **positionOfVirtualParticle;

  int     nX, nY, nZ;

  int     particleNumberOfCenterParticle = 0;
  int     numberOfVirtualParticles;

  double  maxRadius;
  double  minRadius;


  NZERO_setSizeOfVirtualDomain( &nX, &nY, &nZ );

  positionOfVirtualParticle = NZERO_allocateMemoryForVirtualParticles( nX, nY, nZ );

  NZERO_setCoordinateOfVirtualParticles(   positionOfVirtualParticle
										   ,&particleNumberOfCenterParticle
										   ,&numberOfVirtualParticles
										   ,nX, nY, nZ
										   ,averageDistanceBetweenParticles);



  iParticle = particleNumberOfCenterParticle;


  fprintf(FpForLog,"particleNumberOfCenterParticle      = %d\n", particleNumberOfCenterParticle  );
  fprintf(FpForLog,"numberOfVirtualParticles= %d\n", numberOfVirtualParticles);
  fprintf(FpForLog,"nX= %d\n", nX);
  fprintf(FpForLog,"nY= %d\n", nY);
  fprintf(FpForLog,"nZ= %d\n", nZ);


  maxRadius = averageDistanceBetweenParticles * MAXMIN_returnMaxRadius_ratio();
  minRadius = averageDistanceBetweenParticles * MAXMIN_returnMinRadius_ratio();


  NZERO_setZeroToNZeroAndLambda();


  for(jParticle=0; jParticle < numberOfVirtualParticles; jParticle++){

    if(jParticle == iParticle)continue;

    distanceIJ_squared = DISTANCE_calculateSquaredDistance( positionOfVirtualParticle, iParticle, jParticle );
    distanceIJ         = sqrt(distanceIJ_squared);

    if( distanceIJ <= maxRadius ){

      parameter.nZeroOfParticleNumberDensity += WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfParticleNumberDensity );
      /*New*/
      parameter.nZeroOfnumberOfNeighborParticles += WEIGHT_calculateWeightZeroorOneFunction( distanceIJ, parameter.radiusOfParticleNumberDensity );
      parameter.nZeroOfGradient              += WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfGradient );
	  //parameter.nZeroOfGradient              += WEIGHT_calculateWeightFunction( distanceIJ, parameter.maxRadius );
      parameter.nZeroOfLaplacianForPressure  += WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure );
      parameter.nZeroOfLaplacianForViscosity += WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForViscosity);
	  parameter.nZeroOfLaplacianForDiffusion += WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForDiffusion);
	  //parameter.nZeroOfRadius
      parameter.lambdaOfLaplacianForPressure  += distanceIJ_squared * WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure  );
      parameter.lambdaOfLaplacianForViscosity += distanceIJ_squared * WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForViscosity );
      parameter.lambdaTimesNZeroForViscosity  += distanceIJ_squared * WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForViscosity );
	  parameter.lambdaOfLaplacianForDiffusion += distanceIJ_squared * WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForDiffusion);

      if( distanceIJ <= parameter.radiusOfParticleNumberDensity){
	parameter.nZero_withConstantWeight ++;
      }
    }
  }


  parameter.lambdaOfLaplacianForPressure  /= parameter.nZeroOfLaplacianForPressure;
  parameter.lambdaOfLaplacianForViscosity /= parameter.nZeroOfLaplacianForViscosity;
  parameter.lambdaOfLaplacianForDiffusion /= parameter.nZeroOfLaplacianForDiffusion;
  parameter.lambdaOfLaplacianForPressure  *= parameter.relaxationCoefficientOfLambda;

  NZERO_displayCalculatedNZeroAndLambda();

  MEMORY_freeMemoryOf2dimDoubleArray( positionOfVirtualParticle, 3 );

}



void
NZERO_setCoordinateOfVirtualParticles(  
									  double **positionOfVirtualParticle
									  ,int    *particleNumberOfCenterParticle
									  ,int    *numberOfVirtualParticles
									  ,int    nX
									  ,int    nY
									  ,int    nZ
									  ,double averageDistanceBetweenParticles
									  ){
  int     iX, iY, iZ;

  int iParticle = 0;
  (*numberOfVirtualParticles) = 0;

  for(iZ= -nZ; iZ <= nZ; iZ++){
    for(iY= -nY; iY <= nY; iY++){
      for(iX= -nX; iX <= nX; iX++){

		positionOfVirtualParticle[XDIM][iParticle] = iX * averageDistanceBetweenParticles;
		positionOfVirtualParticle[YDIM][iParticle] = iY * averageDistanceBetweenParticles;
		positionOfVirtualParticle[ZDIM][iParticle] = iZ * averageDistanceBetweenParticles;

		if( ((iX==0)&&(iY==0)) &&(iZ==0) ){
		  (*particleNumberOfCenterParticle) = iParticle;
		}
		
		iParticle++;
		(*numberOfVirtualParticles)++;

      }
    }
  }


  if( (*particleNumberOfCenterParticle) == 0){
    OTHER_endProgram("CenterParticle was not able to be detected. [in NZERO_setCoordinateOfVirtualParticles()]");
  }


}



void
NZERO_setZeroToNZeroAndLambda( void ){

  parameter.nZeroOfParticleNumberDensity = 0.0;
  parameter.nZeroOfnumberOfNeighborParticles=0;
  parameter.nZeroOfGradient              = 0.0;
  parameter.nZeroOfLaplacianForPressure  = 0.0; 
  parameter.nZeroOfLaplacianForViscosity = 0.0;
  parameter.nZeroOfLaplacianForDiffusion = 0.0;

  parameter.nZeroOfRigidBody             = 0;

  parameter.lambdaOfLaplacianForPressure  = 0.0;
  parameter.lambdaOfLaplacianForViscosity = 0.0;
  parameter.lambdaTimesNZeroForViscosity  = 0.0;
  parameter.lambdaOfLaplacianForDiffusion = 0.0;

  parameter.nZero_withConstantWeight = 0.0;

}



void
NZERO_setSizeOfVirtualDomain(  
							 int *nX
							 ,int *nY
							 ,int *nZ
							 ){

  double maxRadius_ratio = MAXMIN_returnMaxRadius_ratio();

  (*nX) = (int)(maxRadius_ratio + 2);
  (*nY) = (int)(maxRadius_ratio + 2);
  (*nZ) = (int)(maxRadius_ratio + 2);

  if(NumberOfDimensions == 2){
    (*nZ) = 0;
  }


}


double**
NZERO_allocateMemoryForVirtualParticles( int nX, int nY, int nZ ){

  int    sizeOfArray;
  double **positionOfVirtualParticle;
  
  sizeOfArray = (2 * nX + 1) * (2 * nY + 1) * (2 * nZ + 1);

  positionOfVirtualParticle = MEMORY_allocateMemoryFor2dimDoubleArray( 3, sizeOfArray
							     ,"positionOfVirtualParticle [in NZERO_allocateMemoryForVirtualCoordinate()]");

  return(positionOfVirtualParticle);

}





void
NZERO_displayCalculatedNZeroAndLambda( void ){

  fprintf(FpForLog,"\n\n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        nZero    [in nzero.c]                       \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"parameter.nZeroOfParticleNumberDensity  = %.10e\n", parameter.nZeroOfParticleNumberDensity );
  fprintf(FpForLog,"parameter.nZeroOfGradient               = %.10e\n", parameter.nZeroOfGradient              );
  fprintf(FpForLog,"parameter.nZeroOfLaplacianForPressure   = %.10e\n", parameter.nZeroOfLaplacianForPressure  );
  fprintf(FpForLog,"parameter.nZeroOfLaplacianForViscosity  = %.10e\n", parameter.nZeroOfLaplacianForViscosity );
  fprintf(FpForLog,"parameter.nZeroOfLaplacianForDiffusion  = %.10e\n", parameter.nZeroOfLaplacianForDiffusion);
  fprintf(FpForLog,"parameter.nZero_withConstantWeight      = %.10e\n", parameter.nZero_withConstantWeight );

  fprintf(FpForLog,"parameter.lambdaOfLaplacianForPressure  = %.10e\n", parameter.lambdaOfLaplacianForPressure );
  fprintf(FpForLog,"parameter.lambdaOfLaplacianForViscosity = %.10e\n", parameter.lambdaOfLaplacianForViscosity);
  fprintf(FpForLog,"parameter.lambdaTimesNZeroForViscosity  = %.10e\n", parameter.lambdaTimesNZeroForViscosity );
  fprintf(FpForLog,"parameter.lambdaOfLaplacianForDiffusion = %.10e\n", parameter.lambdaOfLaplacianForDiffusion);
  fprintf(FpForLog,"\n\n");

}




