#include <stdio.h>
#include <math.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "neigh.h"
#include "other.h"
#include "collision.h"

void
COLLISION_calculateCollisionBetweenParticles( void ){


	int    iParticle,jParticle;
	int    iDim;
	int    iNeigh;   

	double vectorIJ[3];

	double distanceIJ;
	double distanceIJ_squared;

	double massOfParticleI;  
	double massOfParticleJ;  

	double normalVectorIJ[3]; 
	double impulse;      

	double velocityIMinusVelocityJ[3]; 
	                 
	double velocityCorrectionOfParticleI[3];
	double velocityCorrectionOfParticleJ[3];

	double positionCorrectionOfParticleI[3];
	double positionCorrectionOfParticleJ[3];
	double dot;

	int countOfCollisionBetween_fluidAndFluid     = 0;
	int countOfCollisionBetween_fluidAndWall      = 0;
	int countOfCollisionBetween_fluidAndDummyWall = 0;
	int countOfCollision                          = 0;

	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	NEIGH_selectNeighborTable(   &numberOfNeighborParticles
				     ,&neighborTable
				     ,parameter.collisionDistance_ratio
					 ,&neighborTablePeriodic
				     );

	/*===================================================================
	  -force * normalVectorIJ * dt = massOfParticleI * velocityI' + massOfParticleI * velocityI
	   force * normalVectorIJ * dt = massOfParticleI * velocityJ' + massOfParticleJ * velocityJ
	   collisionCoefficient = - (velocityI' dot normalVectorIJ - velocityJ' dot normalVectorIJ)
	                        /(velocityI dot normalVectorIJ - velocityJ dot normalVectorIJ)
          ===================================================================*/

#pragma omp parallel for private(jParticle, iDim, iNeigh, vectorIJ, distanceIJ, distanceIJ_squared, massOfParticleI, massOfParticleJ, normalVectorIJ, impulse, velocityIMinusVelocityJ, velocityCorrectionOfParticleI, velocityCorrectionOfParticleJ, positionCorrectionOfParticleI, positionCorrectionOfParticleJ, dot) shared(countOfCollisionBetween_fluidAndFluid, countOfCollisionBetween_fluidAndWall, countOfCollisionBetween_fluidAndDummyWall, countOfCollision)
	for(iParticle = 0;  iParticle < particle.totalNumber;  iParticle++) {

	  if(particle.type[iParticle] == GHOST)continue;
	  massOfParticleI = physicalProperty.massDensity[particle.type[iParticle]];


	  for(iNeigh=0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
	    jParticle = neighborTable[iParticle][iNeigh];

	    if(jParticle <= iParticle)continue;
	    if(particle.type[jParticle]==GHOST)continue;

	    distanceIJ_squared = 0.0;

	    for(iDim=0; iDim < NumberOfDimensions; iDim++){
	      vectorIJ[iDim] = particle.position[iDim][jParticle] - particle.position[iDim][iParticle];
	      if((iDim==0) && (neighborTablePeriodic[iParticle][iNeigh]==1))
			  vectorIJ[iDim] -= parameter.pipeLength;
		  else if((iDim==0) && (neighborTablePeriodic[iParticle][iNeigh]==2))
			  vectorIJ[iDim] += parameter.pipeLength;
	      distanceIJ_squared   += vectorIJ[iDim] * vectorIJ[iDim];
	    }


	    if(distanceIJ_squared < parameter.collisionDistance_squared) {

	      distanceIJ = sqrt(distanceIJ_squared);

	      massOfParticleJ = physicalProperty.massDensity[particle.type[jParticle]];

	      for(iDim = 0; iDim < NumberOfDimensions; iDim++){
		normalVectorIJ[iDim] = vectorIJ[iDim]/distanceIJ;
	      }


	      for(iDim = 0; iDim < NumberOfDimensions; iDim++){
		velocityIMinusVelocityJ[iDim] = particle.velocity[iDim][iParticle] - particle.velocity[iDim][jParticle];
	      }

	      /* ---------------------------------------------------------
		    force * dt = ( 1 + collisionCoefficient ) ( velocityI - velocityJ ) 
                               dot normaoVectorIJ (massOfParticleI * massOfParticleJ)
                               / (massOfParticleI + massOfParticleJ)
                 --------------------------------------------------------*/
	      
	      dot = 0.0;
	      for(iDim = 0; iDim < NumberOfDimensions; iDim++){
		dot += velocityIMinusVelocityJ[iDim] * normalVectorIJ[iDim];
	      }

	      impulse =   ( 1.0 + parameter.collisionCoefficient) * dot * (massOfParticleI * massOfParticleJ) /(massOfParticleI + massOfParticleJ);


	      /* -------------------------------------------------------------- 
		 If the two particles are separating, collision does not occur,
		 even if the distance between particles is too short.
		 --------------------------------------------------------------*/
	      if(impulse < 0.0) continue;

	      for(iDim = 0; iDim < NumberOfDimensions; iDim++){	      
			velocityCorrectionOfParticleI[iDim] = - (impulse / massOfParticleI) * normalVectorIJ[iDim];
			velocityCorrectionOfParticleJ[iDim] = + (impulse / massOfParticleJ) * normalVectorIJ[iDim];

			positionCorrectionOfParticleI[iDim] =   velocityCorrectionOfParticleI[iDim] * timer.dt;
			positionCorrectionOfParticleJ[iDim] =   velocityCorrectionOfParticleJ[iDim] * timer.dt;
	      }

	      /*
	      COLLISION_displayThatCollisionOccured( pcl, iParticle, jParticle, distanceIJ, impulse);
	      */

	      COLLISION_countCollision(  &countOfCollision
					 ,&countOfCollisionBetween_fluidAndFluid
					 ,&countOfCollisionBetween_fluidAndWall
					 ,&countOfCollisionBetween_fluidAndDummyWall
					 ,iParticle
					 ,jParticle
					 );


	      if((particle.type[iParticle] != parameter.dummyWallType)&&(particle.type[iParticle] != parameter.wallType)){
		for(iDim = 0; iDim < NumberOfDimensions; iDim++){	       
		  particle.velocity[iDim][iParticle] += velocityCorrectionOfParticleI[iDim];
		  particle.position[iDim][iParticle] += positionCorrectionOfParticleI[iDim];
		}
	      }


	      if((particle.type[jParticle] != parameter.dummyWallType)&&(particle.type[jParticle] != parameter.wallType)){

		for(iDim = 0; iDim < NumberOfDimensions; iDim++){	       
		  particle.velocity[iDim][jParticle] += velocityCorrectionOfParticleJ[iDim];
		  particle.position[iDim][jParticle] += positionCorrectionOfParticleJ[iDim];
		}
	      }

	    }
	  }
	}

	
	COLLISION_displayCountOfCollision(  countOfCollision
					    ,countOfCollisionBetween_fluidAndFluid
					    ,countOfCollisionBetween_fluidAndWall
					    ,countOfCollisionBetween_fluidAndDummyWall
					    );

}




void
COLLISION_countCollision(
						 int *countOfCollision
						 ,int *countOfCollisionBetween_fluidAndFluid
						 ,int *countOfCollisionBetween_fluidAndWall
						 ,int *countOfCollisionBetween_fluidAndDummyWall
						 ,int  iParticle
						 ,int  jParticle
						 ){

  if(( particle.type[iParticle] != parameter.wallType ) && ( particle.type[iParticle] != parameter.dummyWallType )){

	if(( particle.type[jParticle] != parameter.wallType ) && ( particle.type[jParticle] != parameter.dummyWallType )){

	  (*countOfCollisionBetween_fluidAndFluid)++;

	}else if( particle.type[jParticle] == parameter.wallType ){

	  (*countOfCollisionBetween_fluidAndWall)++;

	}else if( particle.type[jParticle] == parameter.dummyWallType ){

	  (*countOfCollisionBetween_fluidAndDummyWall)++;
	}


  }else{

	if( particle.type[iParticle] == parameter.wallType ){

	  (*countOfCollisionBetween_fluidAndWall)++;

	}else if( particle.type[iParticle] == parameter.dummyWallType ){

	  (*countOfCollisionBetween_fluidAndDummyWall)++;

	}
  }

  (*countOfCollision)++;

}




void
COLLISION_displayCountOfCollision( 
								   int countOfCollision
								  ,int countOfCollisionBetween_fluidAndFluid
								  ,int countOfCollisionBetween_fluidAndWall
								  ,int countOfCollisionBetween_fluidAndDummyWall
								  ){

  if( countOfCollision >= 1){
    fprintf(FpForLog,"Collision occured. -----total:               %d\n", countOfCollision );
    fprintf(FpForLog,"                        fluid and fluid:     %d\n", countOfCollisionBetween_fluidAndFluid     );
    fprintf(FpForLog,"                        fluid and wall:      %d\n", countOfCollisionBetween_fluidAndWall      );
    fprintf(FpForLog,"                        fluid and dummyWall: %d\n", countOfCollisionBetween_fluidAndDummyWall );
  }

}



void
COLLISION_displayThatCollisionOccured(
										 int            iParticle
										,int            jParticle
										,double         distanceIJ
										,double         impulse
										){

  fprintf(FpForLog,"Collision occured between %d-th particle and %d-th particle.\ntype[%d]=%d, type[%d]=%d  DistanceBetweenParticles= %lf[m]  impulse=%lf [Newton * sec]\n"
		      ,iParticle, jParticle, iParticle, particle.type[iParticle],jParticle,particle.type[jParticle], distanceIJ, impulse );

}
