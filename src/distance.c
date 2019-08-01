#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "other.h"
#include "distance.h"




void
DISTANCE_transformUnitOfDistanceFromRatioToMeter( void ){

  parameter.radiusOfParticleNumberDensity = parameter.radiusOfParticleNumberDensity_ratio * particle.averageDistance;

  parameter.radiusOfGradient              = parameter.radiusOfGradient_ratio              * particle.averageDistance;

  parameter.radiusOfGradient_squared      = parameter.radiusOfGradient * parameter.radiusOfGradient;

  parameter.radiusOfLaplacianForViscosity = parameter.radiusOfLaplacianForViscosity_ratio * particle.averageDistance;

  parameter.radiusOfLaplacianForPressure  = parameter.radiusOfLaplacianForPressure_ratio  * particle.averageDistance;

  parameter.radiusOfLaplacianForDiffusion = parameter.radiusOfLaplacianForDiffusion_ratio * particle.averageDistance;

  parameter.collisionDistance         = parameter.collisionDistance_ratio * particle.averageDistance;

  parameter.collisionDistance_squared = parameter.collisionDistance       * parameter.collisionDistance;

}




double
DISTANCE_calculateDistanceBetweenParticles(
										   double **position
										   ,int   iParticle
										   ,int   jParticle
										   ){

  double distanceIJ = 0.0;
  double xij,yij,zij;


  xij = (position[XDIM][jParticle] - position[XDIM][iParticle]);
  yij = (position[YDIM][jParticle] - position[YDIM][iParticle]);

  if(NumberOfDimensions == 3){

    zij = (position[ZDIM][jParticle] - position[ZDIM][iParticle]);
    distanceIJ = sqrt( xij*xij + yij*yij + zij*zij);        

  }else{
    distanceIJ = sqrt( xij*xij + yij*yij);        
  }
   
  return(distanceIJ);

}

double
DISTANCE_calculateDistanceBetweenParticlesPeriodic(
										   double **position
										   ,int   iParticle
										   ,int   jParticle
										   ,int   Boundary
										   ){

  double distanceIJ = 0.0;
  double xij,yij,zij;


  xij = (position[XDIM][jParticle] - position[XDIM][iParticle]);
  if(Boundary==1)
	  xij -= parameter.pipeLength;
  else if(Boundary==2)
	  xij += parameter.pipeLength;
  yij = (position[YDIM][jParticle] - position[YDIM][iParticle]);

  if(NumberOfDimensions == 3){

    zij = (position[ZDIM][jParticle] - position[ZDIM][iParticle]);
    distanceIJ = sqrt( xij*xij + yij*yij + zij*zij);        

  }else{
    distanceIJ = sqrt( xij*xij + yij*yij);        
  }
   
  return(distanceIJ);

}





double
DISTANCE_calculateSquaredDistance(
								  double **position
								  ,int   iParticle
								  ,int   jParticle
								  ){


  double distanceIJ_squared = 0.0;
  double xij,yij,zij;

  xij = (position[XDIM][jParticle] - position[XDIM][iParticle]);
  yij = (position[YDIM][jParticle] - position[YDIM][iParticle]);


  if(NumberOfDimensions == 3){
    zij = (position[ZDIM][jParticle] - position[ZDIM][iParticle]);
    distanceIJ_squared = (xij*xij + yij*yij + zij*zij);        
  }else{
    distanceIJ_squared = (xij*xij + yij*yij);        
  }

  return( distanceIJ_squared );

}

double
DISTANCE_calculateSquaredDistancePeriodic(
								  double **position
								  ,int   iParticle
								  ,int   jParticle
								  ,int   Boundary
								  ){


  double distanceIJ_squared = 0.0;
  double xij,yij,zij;

  xij = (position[XDIM][jParticle] - position[XDIM][iParticle]);
  if(Boundary==1)
	  xij -= parameter.pipeLength;
  else if(Boundary==2)
	  xij += parameter.pipeLength;
  yij = (position[YDIM][jParticle] - position[YDIM][iParticle]);


  if(NumberOfDimensions == 3){
    zij = (position[ZDIM][jParticle] - position[ZDIM][iParticle]);
    distanceIJ_squared = (xij*xij + yij*yij + zij*zij);        
  }else{
    distanceIJ_squared = (xij*xij + yij*yij);        
  }

  return( distanceIJ_squared );

}





double
DISTANCE_calculateDistanceBetweenTwoPositions( double *position1, double *position2 ){

  int    iDim;
  double distance = 0.0;

  for(iDim=0; iDim < NumberOfDimensions; iDim++) {
	distance += pow( position2[iDim] - position1[iDim], 2.0 );
  }

  distance = sqrt(distance);

  return distance;

}




