#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "other.h"
#include "pressure.h"
#include "memory.h"
#include "mathe.h"
#include "copy.h"
#include "other.h"
#include "solver.h"

void
SOLVER_solveSimultaniousEquations( void ){

  double smallNumberForCheckingConvergence;

  smallNumberForCheckingConvergence = SOLVER_returnSmallNumberForCheckingConvergence();

  parameter.relaxationCoefficientForSorSolver = 1.2;

  SOLVER_setInitialSolution();

  SOLVER_conjugateGradientMethod(  particle.totalNumber
				   ,particle.coefficientMatrixOfPressure
				   ,particle.pressure
				   ,particle.sourceTermOfPressure
				   ,particle.flagOfBoundaryCondition
				   ,particle.numberOfNeighborParticles_large
				   ,particle.neighborTable_large
				   ,parameter.maxIterationNumberInIterationSolver
				   ,parameter.minIterationNumberInIterationSolver
				   ,smallNumberForCheckingConvergence
				   );
  omp_set_num_threads(parameter.threads);

}





double
SOLVER_returnSmallNumberForCheckingConvergence( void ){

  double epsilon;

  epsilon = parameter.smallNumberForCheckingConvergenceInIterationSolver;

  /*
  epsilon = parameter.smallNumberForCheckingConvergenceInIterationSolver;
  */
  /*
  epsilon = (timer.dt * timer.dt) * parameter.smallNumberForCheckingConvergenceInIterationSolver;
  */
  /*
  epsilon = (timer.dt * particle.averageDistance) * parameter.smallNumberForCheckingConvergenceInIterationSolver;
  */

  return epsilon;

}


void
SOLVER_setInitialSolution( void ){

  OTHER_fillOneDimDoubleArrayWithZero( particle.totalNumber, particle.pressure );

}



int
SOLVER_conjugateGradientMethod(
			       int      totalNumber
			       ,double **matrix
			       ,double  *solution
			       ,double  *sourceTerm
			       ,int     *flagOfBoundaryCondition
			       ,int     *numberOfNeighborParticles
			       ,int     **neighborTable
			       ,int      maxIteration
			       ,int      minIteration
			       ,double   smallNumberForCheckingConvergence
				){

    double  *solution_next;
    double  *remainder;
    double  *remainder_next;

    double  *direction;
    double  *direction_next;
    double  *matrixTimesDirection;
    double  *matrixTimesSolution;

    double   beta;
    double   alpha;
    double   dot1,dot2;

    int      iIteration;
	double   normOfSourceTerm;

    solution_next  = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber, "solution_next [in solver.c]" );
    remainder      = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber, "remainder     [in solver.c]" );
    remainder_next = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber, "remainder_next[in solver.c]" );
    direction      = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber, "direction     [in solver.c]" );
    direction_next = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber, "direction_next[in solver.c]" );

    matrixTimesDirection = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber
								    ,"matrixTimesDirection[in solver.c]" );

    matrixTimesSolution = MEMORY_allocateMemoryFor1dimDoubleArray( totalNumber
								   ,"matrixTimesSolution[in solver.c]" );

    normOfSourceTerm = SOLVER_calculateNormOfSorceTerm(  totalNumber
							 ,sourceTerm
							 ,flagOfBoundaryCondition
							 );

    OTHER_fillOneDimDoubleArrayWithZero( totalNumber, matrixTimesSolution  );

    OTHER_fillOneDimDoubleArrayWithZero( totalNumber, matrixTimesDirection );

    MATHE_multiplyMatrixByVector( totalNumber, matrixTimesSolution, matrix, solution, numberOfNeighborParticles, neighborTable, flagOfBoundaryCondition);

    MATHE_subtractVectorFromVector( totalNumber, remainder, sourceTerm, matrixTimesSolution, 1.0);

    COPY_copy1dimDoubleArray( totalNumber, direction, remainder );

    for( iIteration = 0; iIteration < maxIteration; iIteration++){

      dot1 = MATHE_innerProductOfVectorsForIterationSolver( remainder, remainder, totalNumber,  flagOfBoundaryCondition);

      MATHE_multiplyMatrixByVector( totalNumber, matrixTimesDirection, matrix, direction, numberOfNeighborParticles, neighborTable, flagOfBoundaryCondition);

      dot2 = MATHE_innerProductOfVectorsForIterationSolver( direction, matrixTimesDirection, totalNumber,  flagOfBoundaryCondition);

      alpha = dot1 / dot2;

      MATHE_addVectorToVector( totalNumber, solution_next, solution,  direction, alpha);

      MATHE_subtractVectorFromVector( totalNumber, remainder_next, remainder,  matrixTimesDirection, alpha);

	  /*
      if ( YES == SOLVER_checkConvergence( totalNumber, flagOfBoundaryCondition, remainder_next, iIteration, minIteration, smallNumberForCheckingConvergence )){
	  */

      if ( YES == SOLVER_checkConvergence( totalNumber, flagOfBoundaryCondition, remainder_next, normOfSourceTerm, iIteration, minIteration, smallNumberForCheckingConvergence )){

	COPY_copy1dimDoubleArray( totalNumber, solution, solution_next );

	SOLVER_displayStateOfConvergence( FpForLog, iIteration, minIteration );

	free(solution_next);
	free(remainder);
	free(remainder_next);
	free(direction);
	free(direction_next);
	free(matrixTimesDirection);
	free(matrixTimesSolution);
        return(YES);

      }


      dot2 = MATHE_innerProductOfVectorsForIterationSolver( remainder_next, remainder_next, totalNumber,  flagOfBoundaryCondition);

      beta = dot2 / dot1;

      MATHE_addVectorToVector( totalNumber, direction_next, remainder_next, direction, beta);

      COPY_copy1dimDoubleArray( totalNumber, remainder, remainder_next);
      COPY_copy1dimDoubleArray( totalNumber, direction, direction_next);
      COPY_copy1dimDoubleArray( totalNumber, solution,  solution_next );

    }

    SOLVER_displayWarnigMessageForExcessOfIterationNumber( iIteration, maxIteration );

    free(solution_next);
    free(remainder);
    free(remainder_next);
    free(direction);
    free(direction_next);
    free(matrixTimesDirection);
    free(matrixTimesSolution);

    return(NO);

}


void
SOLVER_displayStateOfConvergence( FILE *fp, int iIteration, int minIteration ){

  fprintf(fp,"Number of iterations in solver of pressure: %10d\n",iIteration);

  if(iIteration == minIteration){
    fprintf(fp,  "WARNING: Number of iterations is too few.  The calculated pressure might be zero or wrong.\n");
  }
}



void
SOLVER_displayWarnigMessageForExcessOfIterationNumber( int iIteration, int maxIteration ){

  fprintf(FpForLog,"Number of iterations in solver of pressure: %10d\n",iIteration);
  fprintf(FpForLog,"WARNING: Number of iterations exceeded the limit.\n");
  fprintf(FpForLog,"         The upper limit of iteration number = %d\n", maxIteration);


  fprintf(stderr,"Number of iterations in solver of pressure: %10d\n",iIteration);
  fprintf(stderr,"WARNING: Number of iterations exceeded the limit.\n");
  fprintf(stderr,"         The upper limit of iteration number = %d\n", maxIteration);

}



double
SOLVER_calculateNormOfSorceTerm(
				int    totalNumber
				,double *sourceTerm
				,int    *flagOfBoundaryCondition
				  ){

  int    iParticle;
  int    count;
  double normOfSourceTerm;

  count            = 0;
  normOfSourceTerm = 0.0;

#pragma omp parallel for reduction(+:normOfSourceTerm)
  for(iParticle=0; iParticle < totalNumber; iParticle++){

    if( flagOfBoundaryCondition[iParticle] != INNER_PARTICLE   )continue;

    normOfSourceTerm += pow(sourceTerm[iParticle], 2.0);
  }

  return sqrt(normOfSourceTerm);

}


int
SOLVER_checkConvergence(
			int    totalNumber
			,int    *flagOfBoundaryCondition
			,double *remainder
			,double normOfSourceTerm
			,int    iIteration
			,int    minIteration
			,double smallNumberForCheckingConvergence
			  ){

  int    iParticle;
  double sumOfRemainder;
  int    count;

  count = 0;
  sumOfRemainder  = 0.0;

  for(iParticle=0; iParticle < totalNumber; iParticle++){


    sumOfRemainder  += pow(remainder[iParticle], 2.0);
    count++;

  }


  if ( ( sqrt(sumOfRemainder) <=  smallNumberForCheckingConvergence * normOfSourceTerm  ) &&  (iIteration >= minIteration )) {
    return(YES);
  }else{
    return(NO);
  }


}



int
SOLVER_checkConvergence_NORMAL(
			       int    totalNumber
			       ,int    *flagOfBoundaryCondition
			       ,double *remainder
			       ,int    iIteration
			       ,int    minIteration
			       ,double smallNumberForCheckingConvergence
				 ){

  int    iParticle;
  double absoluteValueOfRemainder;
  double sumOfRemainders;
  int    count;

  count = 0;
  sumOfRemainders = 0.0;

  for(iParticle=0; iParticle < totalNumber; iParticle++){

    if( flagOfBoundaryCondition[iParticle] != INNER_PARTICLE   )continue;


    absoluteValueOfRemainder = fabs(remainder[iParticle]);

    if(absoluteValueOfRemainder > smallNumberForCheckingConvergence){
      count++;
    }

    sumOfRemainders += absoluteValueOfRemainder;

  }


  if ( (count==0) &&  (iIteration >= minIteration )) {
    return(YES);
  }else{
    return(NO);
  }


}



int
SOLVER_findIndexOfIParticleInNeighborListOfJParticle( int iParticle,  int jParticle){

  int iNeigh;
  int index = -1;

  for(iNeigh=0; iNeigh < particle.numberOfNeighborParticles_large[jParticle]; iNeigh++){
    if(iParticle == particle.neighborTable_large[jParticle][iNeigh]){
      index = iNeigh;
      return index;
    }
  }

  if( index < 0){
    OTHER_endProgram("ERROR: in findIndexOfIParticleInNeighborListOfJParticle()\n");
  }

  return index;

}


