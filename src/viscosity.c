#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "distance.h"
#include "weight.h"
#include "neigh.h"
#include "other.h"
#include "viscosity.h"

#include "mathe.h"
#include "density.h"
#include "solver.h"
#include "copy.h"


void
VISCOSITY_calculateViscosity( void ){

  int iParticle, jParticle;
  int iDim;
  int iNeigh;   

  double sigma[MAX_DIM];
  double distanceIJ;
  double weightIJ;

  int    *numberOfNeighborParticles;
  int    **neighborTable;
  int    **neighborTablePeriodic;
  double kinematicViscosity;
  double radiusOfLaplacianForViscosity;
  double lambdaTimesNZero;
  double dt;


  if(parameter.flagOfViscosityCalculation == OFF) return;


  kinematicViscosity            = physicalProperty.kinematicViscosity;
  radiusOfLaplacianForViscosity = parameter.radiusOfLaplacianForViscosity;
  lambdaTimesNZero             = parameter.lambdaTimesNZeroForViscosity;
  dt                            = timer.dt;


  NEIGH_selectNeighborTable(  &numberOfNeighborParticles
							  ,&neighborTable
							  ,parameter.radiusOfLaplacianForViscosity_ratio
							  ,&neighborTablePeriodic
							  );


#pragma omp parallel for private(jParticle,iDim,iNeigh,sigma,distanceIJ,weightIJ)
     for( iParticle = 0; iParticle < particle.totalNumber; iParticle++){

          if(particle.type[iParticle]==GHOST                  ) continue;
          if(particle.type[iParticle]==parameter.dummyWallType) continue;
          if(particle.type[iParticle]==parameter.wallType     ) continue;

          for(iDim = 0; iDim < NumberOfDimensions; iDim++){
              sigma[iDim] = 0.0;
          }

          for(iNeigh=0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
              jParticle = neighborTable[iParticle][iNeigh];

            if( particle.type[jParticle] == GHOST )continue;

	            //distanceIJ = DISTANCE_calculateDistanceBetweenParticles( particle.position, iParticle, jParticle );
				distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);
                weightIJ  = WEIGHT_calculateWeightFunction(distanceIJ, radiusOfLaplacianForViscosity);
            
            if(parameter.flagOfHighViscosityCalculation == ON){
	           for(iDim = 0; iDim < NumberOfDimensions; iDim++){
                   sigma[iDim] += ( particle.velocity[iDim][jParticle] - particle.velocity[iDim][iParticle])*weightIJ;
	           }
            }else{
               for(iDim = 0; iDim < NumberOfDimensions; iDim++){
                    sigma[iDim] += ( particle.velocity_previous[iDim][jParticle] - particle.velocity_previous[iDim][iParticle]) *  weightIJ;
               }
            }

	      }
          for(iDim = 0; iDim < NumberOfDimensions; iDim++){
                 sigma[iDim] *= (2.0 * NumberOfDimensions) / lambdaTimesNZero;
              if(parameter.flagOfHighViscosityCalculation == ON){
                  
                  particle.velocity[iDim][iParticle] += 0.50*kinematicViscosity*sigma[iDim]*dt;
                  
                  
              }else{
                  
                 particle.velocity[iDim][iParticle] += kinematicViscosity*sigma[iDim]*dt;
    
              }
              
          }
    }

    if(parameter.flagOfHighViscosityCalculation == ON){
        
#pragma omp parallel for
        for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
            
            particle.flagOfBoundaryCondition[iParticle]=INNER_PARTICLE;
        }
                
        
        for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    
            VISCOSITY_setSourceTerm( iDim );
    
            VISCOSITY_setCoefficientMatrixOfViscosity( iDim );
    
            VISCOSITY_solveSimultaniousEquations(iDim);
            
            VISCOSITY_setTemporaryVelocityToZero(iDim);
        }
    }
}

void
VISCOSITY_setCoefficientMatrixOfViscosity(int iDim){

	int    iParticle, jParticle;
	int    iNeigh;

	double n0     = parameter.nZeroOfLaplacianForViscosity;
	double lambda = parameter.lambdaOfLaplacianForViscosity;

	double kinematicViscosity = physicalProperty.kinematicViscosity;

	double distanceIJ;
	double coefficientOfParticleJ;

	double weightIJ;


	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	double radiusOfLaplacianForViscosity = parameter.radiusOfLaplacianForViscosity;
	double dt = timer.dt;

	NEIGH_selectNeighborTable(  &numberOfNeighborParticles
		,&neighborTable
		,parameter.radiusOfLaplacianForViscosity_ratio
		,&neighborTablePeriodic
		);

#pragma omp parallel for private(jParticle,iNeigh,distanceIJ,coefficientOfParticleJ,weightIJ)
	for( iParticle=0; iParticle < particle.totalNumber; iParticle++){

		particle.coefficientMatrixOfViscosity[iParticle][0]=1.0;

		for(iNeigh=0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];


			//distanceIJ = DISTANCE_calculateDistanceBetweenParticles( particle.position, iParticle, jParticle );
			distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);
			weightIJ  = WEIGHT_calculateWeightFunction(distanceIJ, radiusOfLaplacianForViscosity);

			coefficientOfParticleJ  = (kinematicViscosity*dt*NumberOfDimensions*weightIJ)/(n0*lambda);

			particle.coefficientMatrixOfViscosity[iParticle][iNeigh+1]  = - coefficientOfParticleJ;
			particle.coefficientMatrixOfViscosity[iParticle][0]        +=   coefficientOfParticleJ;

		}      
	}


}

void
VISCOSITY_setSourceTerm( int iDim ){
    
    int iParticle;
    
#pragma omp parallel for
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
            
        particle.sourceTermOfViscosity[iParticle] = particle.velocity[iDim][iParticle];
            
    }
}

void
VISCOSITY_setTemporaryVelocityToZero( int iDim ){
    
    int iParticle;
    
#pragma omp parallel for
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        
        if((particle.type[iParticle]==GHOST) || (particle.type[iParticle] == parameter.dummyWallType)|| (particle.type[iParticle] == parameter.wallType)){
            
            particle.velocity[iDim][iParticle]=particle.velocity_previous[iDim][iParticle];
            
        }
    }    
}

void
VISCOSITY_solveSimultaniousEquations( int iDim ){
    
    
    VISCOSITY_setInitialSolution();
    
    
    SOLVER_conjugateGradientMethod( particle.totalNumber
                                   ,particle.coefficientMatrixOfViscosity
                                   ,particle.velocity_component
                                   ,particle.sourceTermOfViscosity
                                   ,particle.flagOfBoundaryCondition
                                   ,particle.numberOfNeighborParticles_large
                                   ,particle.neighborTable_large
                                   ,parameter.maxIterationNumberInIterationSolver
                                   ,parameter.minIterationNumberInIterationSolver
                                   ,parameter.smallNumberForCheckingConvergenceInIterationSolver
                                   );
    
    COPY_copy2dimDoubleArrayfrom1dimDoubleArray( iDim, particle.totalNumber, particle.velocity, particle.velocity_component);
}


void
VISCOSITY_setInitialSolution(void){
    
    int iParticle;
    
#pragma omp parallel for
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
        particle.velocity_component[iParticle] = 0.05;
    }
}
