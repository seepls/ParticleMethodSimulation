#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "maxmin.h"
#include "memory.h"
#include "other.h"
#include "bucket.h"
#include "domain.h"

void
DOMAIN_initializeDomain( void ){

  MAXMIN_updateMaxMinOfParticleProperties();

  if(domain.flagOfAutoSettingOfDomainSize == ON){
    DOMAIN_setRangeOfDomain();
  }

}



void
DOMAIN_setRangeOfDomain( void ){

  int    iDim;
  double width[3];

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
	width[iDim] = (particle.maxPosition[iDim] - particle.minPosition[iDim]);
  }


  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    domain.upperLimit_previous[iDim] = domain.upperLimit[iDim];
    domain.lowerLimit_previous[iDim] = domain.lowerLimit[iDim];
  }


  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
	domain.upperLimit[iDim] = particle.maxPosition[iDim] 
	  + (domain.upperMarginRatio[iDim] * width[iDim]);

	domain.lowerLimit[iDim] = particle.minPosition[iDim] 
	  - (domain.lowerMarginRatio[iDim] * width[iDim]);

	domain.upperLimit[iDim] += 2*0.5 *particle.averageDistance;

  }

  parameter.pipeLength=domain.upperLimit_previous[0] - domain.lowerLimit_previous[0] - particle.averageDistance + particle.averageDistance;
  printf("AverageDistance: %f\n",particle.averageDistance);

  DOMAIN_displayUpdatedDomainSize();

}



int
DOMAIN_checkWhetherParticleIsInDomain( int iParticle, double **position ){

  int answer = IN_DOMAIN;
  int iDim;

  if ((particle.type[iParticle] == parameter.dummyWallType) || (particle.type[iParticle] == parameter.wallType))
	  return(answer);

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){	
	if( position[iDim][iParticle] >= domain.upperLimit[iDim]) answer = OUT_OF_DOMAIN;
	if( position[iDim][iParticle] <  domain.lowerLimit[iDim]) answer = OUT_OF_DOMAIN;
	if( TRUE == OTHER_isnan(position[iDim][iParticle])       ) answer = EXCEPTION_OCCURED;
	//if( TRUE == isinf(position[iDim][iParticle])             ) answer = EXCEPTION_OCCURED; 


	if(answer == EXCEPTION_OCCURED){
	  fprintf(FpForLog,"\n");
	  fprintf(FpForLog,"ERROR: The particle position is infinity or NaN (Not a Number)\n");
	  fprintf(FpForLog,"       particle ID = %d\n", iParticle);
	  fflush(FpForLog);


	  fprintf(stderr,"\n");
	  fprintf(stderr,"ERROR: The particle position is infinity or NaN (Not a Number)\n");
	  fprintf(stderr,"       particle ID = %d\n", iParticle);
	  fflush(stderr);

	  fprintf(FpForLog,"\n");
	  fflush(FpForLog);
	  OTHER_endProgram("in func of DOMAIN_checkWhetherParticleIsInDomain()");
	}

  }

  return(answer);

}


void
DOMAIN_relocateParticlePeriodic( int iParticle, double **position ){
	if ((particle.type[iParticle] == parameter.dummyWallType) || (particle.type[iParticle] == parameter.wallType))
	  return;
	
  if( particle.position[0][iParticle] >= domain.upperLimit_previous[0]){
	  //printf("xBefore: %f, ",particle.position[0][iParticle]);
	  particle.position[0][iParticle] -=parameter.pipeLength;
	  //printf("axAfter: %f\n",particle.position[0][iParticle]);
  }
  else if( particle.position[0][iParticle] <  domain.lowerLimit_previous[0]){
	  //printf("xBefore: %f, ",particle.position[0][iParticle]);
	  particle.position[0][iParticle] +=parameter.pipeLength;
	  //printf("xAfter: %f\n",particle.position[0][iParticle]);
  }
  
  
}



void
DOMAIN_displayUpdatedDomainSize( void ){

  int iDim;

  fprintf(FpForLog, "====================================================\n");
  fprintf(FpForLog, "        updated domain size                         \n");
  fprintf(FpForLog, "====================================================\n");

  for( iDim=0; iDim < NumberOfDimensions; iDim++){
	fprintf(FpForLog,"domain.upperLimit[%d]  %lf -> %lf [m]\n"
	    ,iDim, domain.upperLimit_previous[iDim], domain.upperLimit[iDim] );
  }

  for( iDim=0; iDim < NumberOfDimensions; iDim++){
	fprintf(FpForLog,"domain.lowerLimit[%d]  %lf -> %lf [m]\n"
	    ,iDim, domain.lowerLimit_previous[iDim], domain.lowerLimit[iDim] );
  }
  

  fprintf(FpForLog,"\n");	
  fflush(FpForLog);

}
