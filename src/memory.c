#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "memory.h"


int*
MEMORY_allocateMemoryFor1dimIntArray( int num, char *variableName ){

  int *array;

  if( num < 1){

    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
    fprintf(FpForLog,"         The assigned memory size was not adequate.\n");

  }

  array = (int*)malloc( num * sizeof(int));

  if( array == NULL){
    fprintf(FpForLog,"Error: Memory allocation for \"%s\" failed\n",variableName);
    exit(1);
  }

  return array;
}



double*
MEMORY_allocateMemoryFor1dimDoubleArray( int num, char *variableName ){

  double *array;

  if( num < 1){
    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
    fprintf(FpForLog,"         The assigned memory size was not adequate.\n");
  }

  array = (double*)malloc( num * sizeof(double));

  if( array == NULL){
    fprintf(FpForLog,"Error: Memory allocation for \"%s\" failed\n",variableName);
    exit(1);
  }

  return array;

}


    

char*
MEMORY_allocateMemoryFor1dimCharArray( int num, char *variableName ){

  char *array;

  if( num < 1){
    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
    fprintf(FpForLog,"         The assigned memory size was not adequate.\n");
  }

  array = (char*)malloc( num * sizeof(char));

  if( array == NULL){
    fprintf(FpForLog,"Error: Memory allocation for \"%s\" failed\n",variableName);
    exit(1);
  }

  return array;
}


int**
MEMORY_allocateMemoryFor2dimIntArray( int nRow, int nColumn, char *variableName ){

  int **array;
  int iRow;

  if( nRow < 1){
    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
    fprintf(FpForLog,"         The assigned memory size is not adequate.\n");
    fprintf(FpForLog,"         Memory size = %d * %d\n", nRow, nColumn);
  }


  array = (int **)malloc( sizeof(int *) * nRow);

  if (array==NULL){
    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n", variableName);
    fprintf(FpForLog,"         The assigned memory size is not adequate.\n");
    fprintf(FpForLog,"         Memory size = %d * %d\n", nRow, nColumn);
    exit(0);
  }


  for (iRow=0 ; iRow < nRow; iRow++) {
    array[iRow] = (int*)malloc( sizeof(int) * nColumn);

    if (array[iRow]==NULL){
      fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n", variableName);
      fprintf(FpForLog,"         The assigned memory size is not adequate.\n");
      fprintf(FpForLog,"         Memory size = %d * %d\n", nRow, nColumn);
      exit(0);
    }
  }

  return array;

}



double**
MEMORY_allocateMemoryFor2dimDoubleArray( int nRow, int nColumn, char *variableName ){

  double **array;
  int iRow;

  if( nRow < 1){
    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
    fprintf(FpForLog,"         The assigned memory size is not adequate.\n");
    fprintf(FpForLog,"         Memory size = %d * %d\n", nRow, nColumn);
  }


  array = (double **)malloc( sizeof(double *) * nRow);

  if (array==NULL){
    fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
    fprintf(FpForLog,"         The assigned memory size is not adequate.\n");
    fprintf(FpForLog,"         Memory size = %d * %d\n", nRow, nColumn);
    exit(0);
  }

  for (iRow=0 ; iRow < nRow; iRow++) {
    array[iRow] = (double*)malloc( sizeof(double) * nColumn);

    if (array[iRow]==NULL){
      fprintf(FpForLog,"Warning: Memory allocation of \"%s\" failed.\n",variableName);
      fprintf(FpForLog,"         The assigned memory size is not adequate.\n");
      fprintf(FpForLog,"         Memory size = %d * %d\n", nRow, nColumn);
      exit(0);
    }
  }

  return array;

}




void
MEMORY_freeMemoryOf2dimIntArray( int **array, int nLine ){

  int      iLine;

  for(iLine=0; iLine < nLine; iLine++){
    free(array[iLine]);
  }

  free(array);

}


void
MEMORY_freeMemoryOf2dimDoubleArray( double **array, int nLine ){

  int      iLine;

  for(iLine=0; iLine < nLine; iLine++){
    free(array[iLine]);
  }

  free(array);

}



void
MEMORY_allocateMemoryForParticleStructure( void ){

  particle.type              = MEMORY_allocateMemoryFor1dimIntArray(    particle.totalNumber_upperLimit,     "particle.type"            );
  particle.type_initial      = MEMORY_allocateMemoryFor1dimIntArray(    particle.totalNumber_upperLimit,     "particle.type_initial"    );
    
    

  particle.moleOfBubbles    = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit,     "particle.moleOfBubbles"    );
  particle.numberOfBubbles    = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit,     "particle.numberOfBubbles"    );
  particle.concentrationOfImpurities    = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit,     "particle.concentrationOfImpurities"    );

  particle.moleOfBubbles_previous    = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit,     "particle.moleOfBubbles_previous"    );
  particle.numberOfBubbles_previous    = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit,     "particle.numberOfBubbles_previous"    );
  particle.concentrationOfImpurities_previous    = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit,     "particle.concentrationOfImpurities_previous"    );

    

  particle.position          = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.position"         );
  particle.position_previous = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.position_previous");

  particle.velocity            = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity"         );
  particle.velocity_previous   = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity_previous");
  particle.velocity_correction = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity_correction");
    
  if(parameter.flagOfOutputOfTorqueFile == ON ){
      particle.position1          = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.position1"         );
      particle.position2          = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.position2"         );
      particle.position3          = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.position"         );
      particle.position4          = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.position"         );
      particle.velocity1            = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity1"         );
      particle.velocity2            = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity2"         );
      particle.velocity3            = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity3"         );
      particle.velocity4            = MEMORY_allocateMemoryFor2dimDoubleArray( 3, particle.totalNumber_upperLimit, "particle.velocity4"         );
      
  }


  if(parameter.flagOfBiCG == ON){

    particle.velocity_component  = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit*2, "particle.velocity_component");

    particle.pressure          = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit * 2, "particle.pressure"         );
    particle.sourceTermOfPressure   = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit * 2, "particle.sourceTermOfPressure");

    particle.sourceTermOfViscosity   = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit * 2, "particle.sourceTermOfViscosity");

  }else{

    particle.velocity_component  = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.velocity_component");

    particle.pressure          = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.pressure"         );
    particle.sourceTermOfPressure   = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.sourceTermOfPressure");

    particle.sourceTermOfViscosity   = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.sourceTermOfViscosity");

    particle.velocity_component  = MEMORY_allocateMemoryFor1dimDoubleArray(    particle.totalNumber_upperLimit, "particle.velocity_component");

  }

  particle.pressure_previous = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.pressure_previous");

  particle.minPressureAmongNeighbors = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.minPressureAmongNeighbors");

  particle.particleNumberDensity          = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.particleNumberDensity");
  particle.particleNumberDensity_previous = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.particleNumberDensity_previous");
  particle.particleNumberDensity_eachParticleType = MEMORY_allocateMemoryFor2dimDoubleArray( particle.totalNumber_upperLimit,  parameter.numberOfParticleTypes, "particleNumberDensity_eachParticleType");

  if(parameter.flagOfBiCG == ON){
    particle.numberOfNeighborParticles_large = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber_upperLimit * 2, "particle.numOfNeighborParticles_large");

    particle.flagOfBoundaryCondition = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber_upperLimit * 2, "flagOfBoundaryCondition");
  }else{
    particle.numberOfNeighborParticles_large = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber_upperLimit, "particle.numOfNeighborParticles_large");

    particle.flagOfBoundaryCondition = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber_upperLimit, "flagOfBoundaryCondition");
  }

  particle.numberOfNeighborParticles_small = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber_upperLimit * 2, "particle.numOfNeighborParticles_small");
  particle.numberOfNeighborParticles_small = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber_upperLimit, "particle.numOfNeighborParticles_small");


  particle.particleNumberDensity_withConstantWeight = MEMORY_allocateMemoryFor1dimDoubleArray( particle.totalNumber_upperLimit, "particle.particleNumberDensity_withConstantWeight");

  if (parameter.flagOfConcentrationCalculation == ON){
	  particle.concentration = MEMORY_allocateMemoryFor2dimDoubleArray(3, particle.totalNumber_upperLimit, "particle.concentration");
	  particle.concentration_previous = MEMORY_allocateMemoryFor2dimDoubleArray(3, particle.totalNumber_upperLimit, "particle.concentration_previous");
  }

  if (parameter.flagOfInterfaceCalculation == ON){
	  particle.interfaceGradient = MEMORY_allocateMemoryFor2dimDoubleArray(3, particle.totalNumber_upperLimit, "particle.interfaceGradient");
	  particle.interfaceArea = MEMORY_allocateMemoryFor1dimDoubleArray(particle.totalNumber_upperLimit, "particle.interfaceArea");
  }

}



void
MEMORY_allocateMemoryForCoefficientMatrixOfPressure( void ){


  if(parameter.flagOfBiCG == ON){
    particle.coefficientMatrixOfPressure = MEMORY_allocateMemoryFor2dimDoubleArray( particle.totalNumber_upperLimit * 2
										    ,parameter.capacityOfNeighborTable_large + 1, "particle.coefficientMatrixOfPressure"   );
      

particle.coefficientMatrixOfViscosity = MEMORY_allocateMemoryFor2dimDoubleArray( particle.totalNumber_upperLimit * 2
                                            ,parameter.capacityOfNeighborTable_large + 1, "particle.coefficientMatrixOfPressure"   );

      
  }else{
    particle.coefficientMatrixOfPressure = MEMORY_allocateMemoryFor2dimDoubleArray( particle.totalNumber_upperLimit
										    ,parameter.capacityOfNeighborTable_large + 1, "particle.coefficientMatrixOfPressure"   );
      
    
    particle.coefficientMatrixOfViscosity = MEMORY_allocateMemoryFor2dimDoubleArray( particle.totalNumber_upperLimit
                                            ,parameter.capacityOfNeighborTable_large + 1, "particle.coefficientMatrixOfPressure"   );

    
  }

  /*--- where, "+1" is necessary because the initial space is used as the particle's own coefficient. ---*/

}

