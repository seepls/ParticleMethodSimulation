#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "other.h"



void
OTHER_endProgram( char *errorMessage ){

  OTHER_displayErrorMessage(errorMessage);
  exit(EXIT_FAILURE);

}



void
OTHER_displayErrorMessage( char *errorMessage){

  fprintf(FpForLog, "ERROR: %s \n",errorMessage);
  fprintf(stderr,   "ERROR: %s \n",errorMessage);

}


void
OTHER_normalEnd( void ){

  fprintf(FpForLog, "Normal End\n");
  fprintf(stderr,   "Normal End\n");

}



void
OTHER_changeParticleTypeIntoGhost( int iParticle ){

  int iDim;


  OTHER_displayThatParticleBecameGhost( iParticle );

  particle.type[iParticle] = GHOST;

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    particle.position[iDim][iParticle]  = 0.0;
    particle.velocity[iDim][iParticle]  = 0.0;
    particle.velocity_correction[iDim][iParticle] = 0.0;
  }

  particle.pressure[iParticle]      = 0.0;
  particle.particleNumberDensity[iParticle] = 0.0;

}




void
OTHER_displayThatParticleBecameGhost( int iParticle ){

/*
  int iDim;
  fprintf(FpForLog, "WARNING: %d-th particle became GHOST.",iParticle);

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    fprintf(FpForLog, "   particle.position[%d][%d] = %lf", iDim, iParticle, particle.position[iDim][iParticle]);
  }

  fprintf(FpForLog,"\n");
*/

}



void
OTHER_checkThatParticlesAreNotAllGhost( void ){

  int iParticle;
  int flag_allGhost = YES;

#pragma omp parallel for
  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
    if( particle.type[iParticle] != GHOST ) flag_allGhost = NO;
  }

  if(flag_allGhost == YES){
    OTHER_endProgram("All particles are GHOST [in main()]\n  in func of OTHER_checkThatParticlesAreAllGhost\n");
  }


}




void
OTHER_fillOneDimDoubleArrayWithZero( int num, double *array ){

  int i;

#pragma omp parallel for
  for(i=0; i < num; i++){
    array[i] = 0.0;
  }

}


void
OTHER_fill2dimDoubleArrayWithZero( int nRow, int nColumn, double **array ){

  int iRow, jColumn;

  for(iRow=0; iRow < nRow; iRow++){
    for(jColumn=0; jColumn < nColumn; jColumn++){
      array[iRow][jColumn] = 0.0;
    }
  }

}

int
OTHER_isnan( double value ){

	if(value != value){
		return YES;
	}else{
		return NO;
	}


}




int
OTHER_checkWhetherTheTwoPositionsAreNearyEqual( double *position1, double *position2 ){

  int answer = YES;

  double distance;

  if( NumberOfDimensions == 2){

	distance = pow(position1[XDIM] - position2[XDIM], 2.0) + pow(position1[YDIM] - position2[YDIM], 2.0);

  }else{

	distance = pow(position1[XDIM] - position2[XDIM], 2.0) + pow(position1[YDIM] - position2[YDIM], 2.0) + pow(position1[ZDIM] - position2[ZDIM], 2.0);

  }


  distance = sqrt(distance);

  if( distance > 0.001 * particle.averageDistance ){

	answer = NO;

  }

  return answer;


}




