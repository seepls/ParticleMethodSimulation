#include <stdio.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "file.h"
#include "other.h"
#include "maxmin.h"


void MAXMIN_setInitialMaxMinValue( int *iParticle ){

  int jParticle;
  int iDim;

  for(jParticle = 0; jParticle < particle.totalNumber; jParticle++){

    if(particle.type[jParticle] == GHOST) continue;

	/*--- for Tank ---*/
/*
	if( parameter.flagOfTankCalculation == ON){
	  if(particle.position[ZDIM][jParticle] <= 0.0) continue;
	}
*/


    for(iDim = 0; iDim < NumberOfDimensions; iDim++){
      particle.minPosition[iDim] = particle.position[iDim][jParticle];
      particle.maxPosition[iDim] = particle.position[iDim][jParticle];

      particle.minVelocity[iDim] = particle.velocity[iDim][jParticle];
      particle.maxVelocity[iDim] = particle.velocity[iDim][jParticle];
    }

    particle.minPressure = particle.pressure[jParticle];    
    particle.maxPressure = particle.pressure[jParticle];    

    particle.minParticleNumberDensity = particle.particleNumberDensity[jParticle];    
    particle.maxParticleNumberDensity = particle.particleNumberDensity[jParticle];    

    particle.standardParticleNumber = jParticle;

    break;

  }

  if(jParticle >= particle.totalNumber){
    OTHER_endProgram("All particles are GHOST type  [in set_initial_value() ]");
  }
  

  (*iParticle) = jParticle;

}


void
MAXMIN_updateMaxMinOfParticleProperties( void ){

  int iDim;
  int iParticle;
  static int flagOfFirstTime = OFF;


  iParticle = 0;

  MAXMIN_setInitialMaxMinValue( &iParticle );

  for(iParticle = iParticle; iParticle < particle.totalNumber; iParticle++) {

    MAXMIN_compareTheParticleWithMaxMin( iParticle );

  }

  if( flagOfFirstTime == OFF ){
    MAXMIN_displayMaxMinValueOfParticleProperty();

    for(iDim = 0; iDim < NumberOfDimensions; iDim++){
      particle.maxPosition_initial[iDim] = particle.maxPosition[iDim];
      particle.minPosition_initial[iDim] = particle.minPosition[iDim];
    }

  }

  flagOfFirstTime = ON;

}




void MAXMIN_compareTheParticleWithMaxMin( int iParticle ){

  int iDim;

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){

    if( particle.minPosition[iDim] > particle.position[iDim][iParticle] ){
        particle.minPosition[iDim] = particle.position[iDim][iParticle];
    }

    if( particle.maxPosition[iDim] < particle.position[iDim][iParticle] ){
        particle.maxPosition[iDim] = particle.position[iDim][iParticle];
    }


    if( particle.minVelocity[iDim] > particle.velocity[iDim][iParticle] ){
        particle.minVelocity[iDim] = particle.velocity[iDim][iParticle];
    }

    if( particle.maxVelocity[iDim] < particle.velocity[iDim][iParticle] ){
        particle.maxVelocity[iDim] = particle.velocity[iDim][iParticle];
    }

  }

  
  if( particle.minPressure > particle.pressure[iParticle]){
      particle.minPressure = particle.pressure[iParticle];
  }


  if( particle.maxPressure < particle.pressure[iParticle]){
      particle.maxPressure = particle.pressure[iParticle];
  }


  if( particle.minParticleNumberDensity > particle.particleNumberDensity[iParticle]){
      particle.minParticleNumberDensity = particle.particleNumberDensity[iParticle];
  }

  if( particle.maxParticleNumberDensity < particle.particleNumberDensity[iParticle]){
      particle.maxParticleNumberDensity = particle.particleNumberDensity[iParticle];
      particle.standardParticleNumber  = iParticle;
  }

}




void
MAXMIN_calculateMaxAbsoluteValueOfVelocity( void ){

  int    iParticle;
  int    iDim;
  double squaredVelocity;

  particle.maxAbsoluteValueOfVelocity = 0.0;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == GHOST )continue;

    squaredVelocity = 0.0;

    for( iDim = 0; iDim < NumberOfDimensions; iDim++){
      squaredVelocity += particle.velocity[iDim][iParticle] * particle.velocity[iDim][iParticle];
    }

    if( particle.maxAbsoluteValueOfVelocity < squaredVelocity ){
        particle.maxAbsoluteValueOfVelocity = squaredVelocity;
    }

  }

  particle.maxAbsoluteValueOfVelocity = sqrt(particle.maxAbsoluteValueOfVelocity);

  fprintf(FpForLog, "max|velocity|= %lf[m/s]\n", particle.maxAbsoluteValueOfVelocity);

}



void
MAXMIN_setMaxMinOfRadii( void ){

  parameter.maxRadius_ratio   = MAXMIN_returnMaxRadius_ratio();
  parameter.maxRadius         = MAXMIN_returnMaxRadius();
  parameter.maxRadius_squared = parameter.maxRadius * parameter.maxRadius;

  parameter.minRadius_ratio   = MAXMIN_returnMinRadius_ratio();
  parameter.minRadius         = MAXMIN_returnMinRadius();
  parameter.minRadius_squared = parameter.minRadius * parameter.minRadius;

}




double
MAXMIN_returnMaxRadius_ratio( void ){

  double maxRadius_ratio;

  maxRadius_ratio = parameter.radiusOfParticleNumberDensity_ratio;

  if( maxRadius_ratio < parameter.radiusOfGradient_ratio ){
      maxRadius_ratio = parameter.radiusOfGradient_ratio;
  }


  if( maxRadius_ratio < parameter.radiusOfLaplacianForPressure_ratio ){
      maxRadius_ratio = parameter.radiusOfLaplacianForPressure_ratio;
  }

  if( maxRadius_ratio < parameter.radiusOfLaplacianForViscosity_ratio ){
      maxRadius_ratio = parameter.radiusOfLaplacianForViscosity_ratio;
  }

  if (maxRadius_ratio < parameter.radiusOfLaplacianForDiffusion_ratio){
	  maxRadius_ratio = parameter.radiusOfLaplacianForDiffusion_ratio;
  }


  return(maxRadius_ratio);

}



double
MAXMIN_returnMinRadius_ratio( void ){

  double minRadius_ratio;

  minRadius_ratio = parameter.radiusOfParticleNumberDensity_ratio;

  if( minRadius_ratio > parameter.radiusOfGradient_ratio ){
      minRadius_ratio = parameter.radiusOfGradient_ratio;
  }


  if( minRadius_ratio > parameter.radiusOfLaplacianForPressure_ratio ){
      minRadius_ratio = parameter.radiusOfLaplacianForPressure_ratio;
  }

  if( minRadius_ratio > parameter.radiusOfLaplacianForViscosity_ratio ){
      minRadius_ratio = parameter.radiusOfLaplacianForViscosity_ratio;
  }

  if (minRadius_ratio > parameter.radiusOfLaplacianForDiffusion_ratio){
	  minRadius_ratio = parameter.radiusOfLaplacianForDiffusion_ratio;
  }


  return(minRadius_ratio);

}



double
MAXMIN_returnMaxRadius( void ){

  double maxRadius;

  maxRadius  = MAXMIN_returnMaxRadius_ratio();
  maxRadius *= particle.averageDistance;

  return(maxRadius);

}


double
MAXMIN_returnMinRadius( void ){

  double minRadius;

  minRadius  = MAXMIN_returnMinRadius_ratio();
  minRadius *= particle.averageDistance;

  return(minRadius);

}




double
MAXMIN_returnSmallerValue( double value1, double value2){

  double smallerValue;

  smallerValue = value1;

  if( smallerValue > value2){
    smallerValue = value2;
  }

  return(smallerValue);

}



void
MAXMIN_displayMaxMinValueOfParticleProperty( void ){


    int iDim;

    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "=======================================================\n");
    fprintf(FpForLog, "        Max and Min value of initial particle property \n");
    fprintf(FpForLog, "=======================================================\n");


    for(iDim = 0; iDim < NumberOfDimensions; iDim++){
      fprintf(FpForLog, "particle.maxPosition[%s]        %lf [m]\n",FILE_returnDim(iDim), particle.maxPosition[iDim]);
    }
    fprintf(FpForLog, "\n");


    for(iDim = 0; iDim < NumberOfDimensions; iDim++){
      fprintf(FpForLog, "particle.minPosition[%s]        %lf [m]\n",FILE_returnDim(iDim), particle.minPosition[iDim]);
    }
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");


    for(iDim = 0; iDim < NumberOfDimensions; iDim++){
      fprintf(FpForLog, "particle.maxVelocity[%s]        %lf [m/s]\n",FILE_returnDim(iDim), particle.maxVelocity[iDim]);
    }
    fprintf(FpForLog, "\n");


    for(iDim = 0; iDim < NumberOfDimensions; iDim++){
      fprintf(FpForLog, "particle.minVelocity[%s]        %lf [m/s]\n",FILE_returnDim(iDim), particle.minVelocity[iDim]);
    }
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");

    fprintf(FpForLog, "particle.maxPressure                 %lf [Pascal]\n",particle.maxPressure);
    fprintf(FpForLog, "particle.minPressure                 %lf [Pascal]\n",particle.minPressure);
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");

    fprintf(FpForLog, "particle.maxParticleNumberDeinsity           %lf [Particles]\n",particle.maxParticleNumberDensity);
    fprintf(FpForLog, "particle.minParticleNumberDeinsity           %lf [Particles]\n",particle.minParticleNumberDensity);

    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");
    fprintf(FpForLog, "\n");

}

