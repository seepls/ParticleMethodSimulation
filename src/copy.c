#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "copy.h"


void
COPY_setInitialCoodinatesAndVelocities( void ){

  COPY_copy2dimDoubleArray( 3, particle.totalNumber, particle.position_previous, particle.position);
  COPY_copy2dimDoubleArray( 3, particle.totalNumber, particle.velocity_previous, particle.velocity);

}


void
COPY_storeInitialParticleProperty( void ){

  COPY_copy1dimIntArray( particle.totalNumber, particle.type_initial,  particle.type );

}


void
COPY_updateParticleProperty( void ){

  COPY_copy2dimDoubleArray( 3, particle.totalNumber, particle.position_previous,      particle.position);
  COPY_copy2dimDoubleArray( 3, particle.totalNumber, particle.velocity_previous,      particle.velocity);
  COPY_copy1dimDoubleArray( particle.totalNumber, particle.pressure_previous,              particle.pressure);
  COPY_copy1dimDoubleArray( particle.totalNumber, particle.particleNumberDensity_previous, particle.particleNumberDensity);

  COPY_copy1dimDoubleArray( particle.totalNumber, particle.moleOfBubbles_previous, particle.moleOfBubbles);
  COPY_copy1dimDoubleArray( particle.totalNumber, particle.numberOfBubbles_previous, particle.numberOfBubbles);
  COPY_copy1dimDoubleArray( particle.totalNumber, particle.concentrationOfImpurities_previous, particle.concentrationOfImpurities);
  timer.iTimeStep_copy=timer.iTimeStep+1;
}


void
COPY_copy2dimDoubleArray( int nRow, int nColumn, double **overwrittenArray, double **storedArray){
  int iRow, iColumn;

  for( iRow = 0; iRow < nRow; iRow++){
    for(iColumn = 0; iColumn < nColumn; iColumn++){
      overwrittenArray[iRow][iColumn] = storedArray[iRow][iColumn];
    }
  }

}


void
COPY_copy1dimDoubleArray( int num, double *overwrittenArray, double *storedArray ){

  int i;

#pragma omp parallel for
  for( i=0; i < num; i++){
    overwrittenArray[i] = storedArray[i];
  }

}


void
COPY_copy1dimIntArray( int num, int *overwrittenArray, int *storedArray ){

  int i;

  for( i=0; i < num; i++){
    overwrittenArray[i] = storedArray[i];
  }

}

void
COPY_copy1dimDoubleArrayfrom2dimDoubleArray( int nRow, int nColumn, double *overwrittenArray, double **storedArray){
    int iColumn;

        for(iColumn = 0; iColumn < nColumn; iColumn++){
            overwrittenArray[iColumn] = storedArray[nRow][iColumn];
        }

}

void
COPY_copy2dimDoubleArrayfrom1dimDoubleArray( int nRow, int nColumn, double **overwrittenArray, double *storedArray){
    int iColumn;

#pragma omp parallel for
    for(iColumn = 0; iColumn < nColumn; iColumn++){
        overwrittenArray[nRow][iColumn] = storedArray[iColumn];
    }

}

