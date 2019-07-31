#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


#include "define.h"
#include "struct.h"
#include "extern.h"
#include "quaternion.h"


double
MATHE_average( double a1, double a2 ){

  double average;

  average = ( a1 + a2 )/ 2.0;

  return(average);

}


void
MATHE_multiplyMatrixByVector(
			     int    totalNum
			     ,double *answer
			     ,double **matrix
			     ,double *vector
			     ,int    *numberOfNeighborParticles
			     ,int    **neighborList
			     ,int    *flagOfBoundaryCondition
			       ){

	int iParticle,jParticle,iNeigh;

#pragma omp parallel for private(jParticle,iNeigh)
	for(iParticle=0; iParticle < totalNum; iParticle++){

		if ( flagOfBoundaryCondition[iParticle] != INNER_PARTICLE   ) continue;

		if(parameter.flagOfBiCG == ON){

			if( iParticle < particle.totalNumber){

				answer[iParticle] = matrix[iParticle][0] * vector[iParticle + particle.totalNumber];

			}else{

				answer[iParticle] = matrix[iParticle][0] * vector[iParticle - particle.totalNumber];

			}

		}else{

			answer[iParticle] = matrix[iParticle][0] * vector[iParticle];

		}

		/*
		answer[iParticle] = matrix[iParticle][0] * vector[iParticle];
		*/

		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){

			jParticle = neighborList[iParticle][iNeigh];

			if(flagOfBoundaryCondition[jParticle]==GHOST_OR_DUMMY)continue;

			answer[iParticle] += matrix[iParticle][iNeigh+1] * vector[jParticle];

		}
	}

}



void
MATHE_multiplyMatrixByMatrix_matrixSizeIsThree( double answer[3][3]
						,double matrix_left[3][3]
						,double matrix_right[3][3]
						){
  double matrix_tmp[3][3];

  int iRow, iColumn;
  int k;

  for(iRow=0; iRow < 3; iRow++){
	for(iColumn=0; iColumn < 3; iColumn++){

	  matrix_tmp[iRow][iColumn] = 0.0;

	  for(k=0; k < 3; k++){
		matrix_tmp[iRow][iColumn] += matrix_left[iRow][k] * matrix_right[k][iColumn];
	  }
    }
  }


  for(iRow=0; iRow < 3; iRow++){
	for(iColumn=0; iColumn < 3; iColumn++){
	  answer[iRow][iColumn] = matrix_tmp[iRow][iColumn];
	}
  }

}






void
MATHE_addVectorToVector( int num, double *answer, double *vector1, double *vector2, double a){

  int i;

#pragma omp parallel for
  for(i=0; i < num; i++){
    answer[i] = vector1[i] + a * vector2[i];
  }

}


void
MATHE_add2dimArrayTo2dimArray( int nRow, int nColumn, double **answer, double **array1, double **array2 ){

  int iRow, jColumn;

  for(iRow=0; iRow < nRow; iRow++){
    for(jColumn=0; jColumn < nColumn; jColumn++){
      answer[iRow][jColumn] = array1[iRow][jColumn] + array2[iRow][jColumn];
    }
  }


}


void
MATHE_subtractVectorFromVector( int num, double *answer, double *vector1, double *vector2, double a){

  int i;

#pragma omp parallel for
  for(i=0; i < num; i++){
    answer[i] = vector1[i] - a * vector2[i];
  }

}



double
MATHE_innerProductOfVectors( double *vector1, double *vector2, int num){

  double answer;
  int iParticle;

  answer = 0.0;

  for(iParticle=0; iParticle < num; iParticle++){
    answer += (vector1[iParticle] * vector2[iParticle]);
  }

  return(answer);

}




double
MATHE_innerProductOfVectorsForIterationSolver( double *vector1, double *vector2, int num,  int *flagOfBoundaryCondition){

  double answer;
  int iParticle;

  answer = 0.0;

#pragma omp parallel for reduction(+:answer)
  for(iParticle=0; iParticle < num; iParticle++){

    if(flagOfBoundaryCondition[iParticle] == INNER_PARTICLE ){
      answer += (vector1[iParticle] * vector2[iParticle]);
    }

  }

  return(answer);

}



double
MATHE_linearInterpolation(  double value2, double value1,  double m2, double m1, double n2, double n1 ){

  double interpolatedValue;

  double ratio;

  if( m2 == m1){

    interpolatedValue = value1;

  }else if( n2 == n1 ){

    interpolatedValue = 0.5 * (value2 + value1);

  }else{

    ratio = (m2 - m1)/(n2 - n1);

	if( fabs(ratio) < TOO_SMALL_VALUE ){
	  interpolatedValue = value1;
	}else{
	  interpolatedValue = ratio * (value2 - value1) + value1;
	}

  }

  return interpolatedValue;

}



double
MATHE_absoluteValueOfVector( double *vector, int nDimension ){

  int    iDim;
  double absoluteValue;

  absoluteValue = 0.0;

  for(iDim=0; iDim < nDimension; iDim++){
    absoluteValue += (vector[iDim] * vector[iDim]);
  }

  absoluteValue = sqrt(absoluteValue);

  return absoluteValue;

}


void
MATHE_outerProduct(
				    double *answer
				   ,double *vector1
				   ,double *vector2
				   ){

  double answer_tmp[3];

  answer_tmp[XDIM] = (vector1[YDIM] * vector2[ZDIM] - vector1[ZDIM] * vector2[YDIM]);
  answer_tmp[YDIM] = (vector1[ZDIM] * vector2[XDIM] - vector1[XDIM] * vector2[ZDIM]);
  answer_tmp[ZDIM] = (vector1[XDIM] * vector2[YDIM] - vector1[YDIM] * vector2[XDIM]);


  answer[XDIM] = answer_tmp[XDIM];
  answer[YDIM] = answer_tmp[YDIM];
  answer[ZDIM] = answer_tmp[ZDIM];

}


void
MATHE_normalizeVector( double *vector, int nDimension ){

  double absoluteValueOfVector;
  int    iDim;

  absoluteValueOfVector = MATHE_absoluteValueOfVector( vector, nDimension );

  if(absoluteValueOfVector > 0.0){

    for(iDim=0; iDim < nDimension; iDim++){
      vector[iDim] = vector[iDim]/ absoluteValueOfVector;
    }

  }else{

    for(iDim=0; iDim < nDimension; iDim++){
      vector[iDim] = 0.0;
    }
  }

}




void
MATHE_linearTransform( double *answer, double matrix[3][3], double *vector ){

  int iRow, iColumn;
  double answer_tmp[3];


  for( iRow=0; iRow < 3; iRow++ ){
	answer_tmp[iRow] = 0.0;

	for( iColumn=0; iColumn < 3; iColumn++ ){
	  answer_tmp[iRow] += matrix[iRow][iColumn] * vector[iColumn];
	}
  }

  for( iRow=0; iRow < 3; iRow++ ){
	answer[iRow] = answer_tmp[iRow];
  }

}


void
MATHE_matrixInverse2x2(double **matrix, double **inverse){

	double determinant;

	determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

	//printf("det: %f\n", determinant);

	if (fabs(determinant) < TOO_SMALL_VALUE){
		printf("Warning: Determinant is 0 in MATHE_matrixInverse2x2, setting identity matrix");
		inverse[0][0] = 1;
		inverse[0][1] = 0;
		inverse[1][0] = 0;
		inverse[1][1] = 1;
		return;
	}


	inverse[0][0] = matrix[1][1] / determinant;
	inverse[0][1] = - matrix[0][1] / determinant;
	inverse[1][0] = - matrix[1][0] / determinant;
	inverse[1][1] = matrix[0][0] / determinant;

}

