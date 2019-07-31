#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "domain.h"
#include "memory.h"
#include "sort.h"
#include "bucket.h"
#include "maxmin.h"
#include "other.h"
#include "neigh.h"


void
NEIGH_initializeNeighborTable( void ){

  parameter.maxArraySizeOfOccupiedNeighborTable_large = 0;

  parameter.maxArraySizeOfOccupiedNeighborTable_small = 0;

  NEIGH_setCapacityOfNeighborTableAutomatically();

  NEIGH_displayCapacityOfNeighborTable();

  NEIGH_allocateMemoryForNeighborTable();

}




void
NEIGH_setCapacityOfNeighborTableAutomatically( void ){


  double estimatedCapacity_large;

  double estimatedCapacity_small;

  if(parameter.flagOfAutoSettingOfCapacityOfNeighborTable == ON){

    if( NumberOfDimensions == 2 ){

      estimatedCapacity_large = (parameter.maxRadius_ratio + 1.0) * parameter.maxRadius_ratio * M_PI;

      estimatedCapacity_small = (parameter.minRadius_ratio + 1.0) * parameter.minRadius_ratio * M_PI;

    }else{

      estimatedCapacity_large = (4.0/3.0) * M_PI * parameter.maxRadius_ratio * parameter.maxRadius_ratio * parameter.maxRadius_ratio;

      estimatedCapacity_small = (4.0/3.0) * M_PI * parameter.minRadius_ratio * parameter.minRadius_ratio * parameter.minRadius_ratio;

    }

    parameter.capacityOfNeighborTable_large = (int)(estimatedCapacity_large * parameter.marginRatioForSettingCapacityOfLargeNeighborTable);

    parameter.capacityOfNeighborTable_small = (int)(estimatedCapacity_small * parameter.marginRatioForSettingCapacityOfSmallNeighborTable);

    fprintf(FpForLog,"Capacity of large neighborTable= %d\n",parameter.capacityOfNeighborTable_large);
    fprintf(FpForLog,"Capacity of small neighborTable= %d\n",parameter.capacityOfNeighborTable_small);
  }


}




void
NEIGH_setNeighborTable( double **position ){

  int    iParticle, jParticle;
  int    iStoredParticle;
  double xij,yij,zij;

  int iX, iY, iZ;
  int jX, jY, jZ;

  double distanceIJ_squared;

  if (parameter.flagOfPeriodicX == ON){
	  NEIGH_setNeighborTablePeriodic(position);
	  return;
  }

  BUCKET_storeParticlesInBuckets( position );

  NEIGH_resetNeighborTable();

#pragma omp parallel for private(jParticle,iX,iY,iZ,jX,jY,jZ,iStoredParticle,distanceIJ_squared,xij,yij,zij)
  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

	if(particle.type[iParticle] == GHOST)continue;

	BUCKET_findBucketWhereParticleIsStored( &iX, &iY, &iZ, iParticle, position );

	for( jX = iX-1; jX <= iX+1; jX++){
	  if( jX <  0                        ) continue;
	  if( jX >= domain.numberOfBuckets[XDIM]) continue;

	  for( jY = iY-1; jY <= iY + 1; jY++){
		if( jY <  0                        ) continue;
		if( jY >= domain.numberOfBuckets[YDIM]) continue;

		for( jZ = iZ-1; jZ <= iZ + 1; jZ++){
		  if( jZ <  0                        ) continue;
		  if( jZ >= domain.numberOfBuckets[ZDIM]) continue;


		  for( iStoredParticle= 0; iStoredParticle  < domain.bucket[jX][jY][jZ].count; iStoredParticle++){
			jParticle = domain.bucket[jX][jY][jZ].list[iStoredParticle];

			if(jParticle == iParticle)continue;
			if( particle.type[jParticle] == GHOST )continue;
			if((particle.type[jParticle] == parameter.dummyWallType) && (particle.type[iParticle]== parameter.dummyWallType ))continue;

			xij = (position[XDIM][jParticle] - position[XDIM][iParticle]);
			distanceIJ_squared = xij * xij;
			if( distanceIJ_squared > parameter.maxRadius_squared )continue;


			yij = (position[YDIM][jParticle] - position[YDIM][iParticle]);
			distanceIJ_squared += yij * yij;
			if( distanceIJ_squared > parameter.maxRadius_squared )continue;


			if(NumberOfDimensions == 3){
			  zij = (position[ZDIM][jParticle] - position[ZDIM][iParticle]);
			  distanceIJ_squared += zij * zij;
			  if( distanceIJ_squared > parameter.maxRadius_squared )continue;
			}


			NEIGH_checkCapacityOfNeiborTable(iParticle
							 ,particle.numberOfNeighborParticles_large
							 ,particle.neighborTable_large
							 ,parameter.capacityOfNeighborTable_large
							 ,"particle.neighborTable_large");


			NEIGH_addParticleInNeighborTable( iParticle, jParticle
							  ,particle.numberOfNeighborParticles_large
							  ,particle.neighborTable_large );


			if( distanceIJ_squared > parameter.minRadius_squared )continue;


			if( parameter.numberOfNeighborTables > 1){

			  NEIGH_checkCapacityOfNeiborTable(iParticle
							   ,particle.numberOfNeighborParticles_small
							   ,particle.neighborTable_small
							   ,parameter.capacityOfNeighborTable_small
							   ,"particle.neighborTable_small");


			  NEIGH_addParticleInNeighborTable( iParticle, jParticle
							    ,particle.numberOfNeighborParticles_small
							    ,particle.neighborTable_small );

			}
		  }
		}
	  }
	}

  }

  SORT_sortNeighborTable(  particle.totalNumber, particle.numberOfNeighborParticles_large
						  ,particle.neighborTable_large, ASCENDING_ORDER);

  if( parameter.numberOfNeighborTables > 1){
	SORT_sortNeighborTable( particle.totalNumber, particle.numberOfNeighborParticles_small
						   ,particle.neighborTable_small, ASCENDING_ORDER);
  }

  NEIGH_recordMaxArraySizeOfOccupiedNeighborTable();


}


void
NEIGH_setNeighborTablePeriodic( double **position ){

  int    iParticle, jParticle;
  double xij,yij,zij;

  int periodic;

  double distanceIJ_squared;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++)
	  DOMAIN_relocateParticlePeriodic( iParticle, position );

  BUCKET_storeParticlesInBuckets( position );

  NEIGH_resetNeighborTable();

#pragma omp parallel for private(jParticle,xij,yij,zij,periodic,distanceIJ_squared)
  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

	if(particle.type[iParticle] == GHOST)continue;

	for(jParticle=0; jParticle < particle.totalNumber; jParticle++){

			if(jParticle == iParticle)continue;
			if( particle.type[jParticle] == GHOST )continue;
			if((particle.type[jParticle] == parameter.dummyWallType) && (particle.type[iParticle]== parameter.dummyWallType ))continue;

			xij = (position[XDIM][jParticle] - position[XDIM][iParticle]);
			distanceIJ_squared = xij * xij;
			periodic = 0;
			if( distanceIJ_squared > parameter.maxRadius_squared ){
				if( xij > parameter.maxRadius ){xij -= parameter.pipeLength;distanceIJ_squared = xij * xij;periodic = 1;}
				else if( xij < -parameter.maxRadius ){xij += parameter.pipeLength;distanceIJ_squared = xij * xij;periodic = 2;}
			}
			if( distanceIJ_squared > parameter.maxRadius_squared )continue;


			yij = (position[YDIM][jParticle] - position[YDIM][iParticle]);
			distanceIJ_squared += yij * yij;
			if( distanceIJ_squared > parameter.maxRadius_squared )continue;


			if(NumberOfDimensions == 3){
			  zij = (position[ZDIM][jParticle] - position[ZDIM][iParticle]);
			  distanceIJ_squared += zij * zij;
			  if( distanceIJ_squared > parameter.maxRadius_squared )continue;
			}


			NEIGH_checkCapacityOfNeiborTable(iParticle
							 ,particle.numberOfNeighborParticles_large
							 ,particle.neighborTable_large
							 ,parameter.capacityOfNeighborTable_large
							 ,"particle.neighborTable_large");


			NEIGH_addParticleInNeighborTable( iParticle, jParticle
							  ,particle.numberOfNeighborParticles_large
							  ,particle.neighborTable_large );

			particle.neighborTable_largeBoundary[iParticle][particle.numberOfNeighborParticles_large[iParticle] - 1] = periodic;

			if( distanceIJ_squared > parameter.minRadius_squared )continue;


			if( parameter.numberOfNeighborTables > 1){

			  NEIGH_checkCapacityOfNeiborTable(iParticle
							   ,particle.numberOfNeighborParticles_small
							   ,particle.neighborTable_small
							   ,parameter.capacityOfNeighborTable_small
							   ,"particle.neighborTable_small");


			  NEIGH_addParticleInNeighborTable( iParticle, jParticle
							    ,particle.numberOfNeighborParticles_small
							    ,particle.neighborTable_small );

			  particle.neighborTable_smallBoundary[iParticle][particle.numberOfNeighborParticles_small[iParticle] - 1] = periodic;

			}
		  }
  }

  SORT_sortNeighborTable(  particle.totalNumber, particle.numberOfNeighborParticles_large
						  ,particle.neighborTable_large, ASCENDING_ORDER);

  if( parameter.numberOfNeighborTables > 1){
	SORT_sortNeighborTable( particle.totalNumber, particle.numberOfNeighborParticles_small
						   ,particle.neighborTable_small, ASCENDING_ORDER);
  }

  NEIGH_recordMaxArraySizeOfOccupiedNeighborTable();


}




void
NEIGH_recordMaxArraySizeOfOccupiedNeighborTable( void ){

  int iParticle;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == GHOST               )continue;
    if( particle.type[iParticle] == parameter.dummyWallType )continue;

    if( parameter.maxArraySizeOfOccupiedNeighborTable_large < particle.numberOfNeighborParticles_large[iParticle] ){
      parameter.maxArraySizeOfOccupiedNeighborTable_large = particle.numberOfNeighborParticles_large[iParticle];
    }

    if( parameter.numberOfNeighborTables > 1){

      if( parameter.maxArraySizeOfOccupiedNeighborTable_small < particle.numberOfNeighborParticles_small[iParticle] ){
		parameter.maxArraySizeOfOccupiedNeighborTable_small = particle.numberOfNeighborParticles_small[iParticle];
      }

    }

  }

}




void
NEIGH_checkCapacityOfNeiborTable(
								 int  iParticle
								 ,int  *numberOfNeighborParticles
								 ,int  **neighborTable
								 ,int  capacity
								 ,char *nameOfTable
								 ){

  if( numberOfNeighborParticles[iParticle] >= capacity){
    NEIGH_displayErrorMessageForMemoryOverFlow(iParticle, numberOfNeighborParticles, neighborTable, nameOfTable);
    OTHER_endProgram("Neighbor table overflowed. [in NEIGH_checkCapacityOfNeiborTable()]\n Please check the numberOfDimensions and averageDistanceBetweenParticles in dataFile.");
  }

}



void
NEIGH_addParticleInNeighborTable(
								 int iParticle
								 ,int jParticle
								 ,int *numberOfNeighborParticles
								 ,int **neighborTable
								 ){

  int iNeigh;

  iNeigh = numberOfNeighborParticles[iParticle];

  neighborTable[iParticle][iNeigh] = jParticle;
  numberOfNeighborParticles[iParticle]++;

}




void
NEIGH_resetNeighborTable( void ){

  int iParticle;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
    particle.numberOfNeighborParticles_large[iParticle] = 0;
  }


  if( parameter.numberOfNeighborTables > 1 ){

	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
	  particle.numberOfNeighborParticles_small[iParticle] = 0;
	}

  }

}




void
NEIGH_selectNeighborTable(  int **numberOfNeighborParticles
							,int ***neighborTable
							,double radius_ratio
							,int ***neighborTablePeriodic
							){


  double minRadius_ratio;

  minRadius_ratio   = MAXMIN_returnMinRadius_ratio();

  if( parameter.numberOfNeighborTables == 1 ){

	(*numberOfNeighborParticles) = particle.numberOfNeighborParticles_large;
	(*neighborTable)             = particle.neighborTable_large;
	(*neighborTablePeriodic)     = particle.neighborTable_largeBoundary;

  }else if( radius_ratio > minRadius_ratio ){

	(*numberOfNeighborParticles) = particle.numberOfNeighborParticles_large;
	(*neighborTable)             = particle.neighborTable_large;
	(*neighborTablePeriodic)     = particle.neighborTable_largeBoundary;

  }else{

	(*numberOfNeighborParticles) = particle.numberOfNeighborParticles_small;
	(*neighborTable)             = particle.neighborTable_small;
	(*neighborTablePeriodic)     = particle.neighborTable_smallBoundary;

  }

}




void
NEIGH_displayCapacityOfNeighborTable( void ){

  fprintf(FpForLog,"\n\n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        Capacity of neighbor table                  \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"parameter.capacityOfNeighborTable_large = %d\n", parameter.capacityOfNeighborTable_large );

  if( parameter.numberOfNeighborTables > 1 ){
	fprintf(FpForLog,"parameter.capacityOfNeighborTable_small = %d\n", parameter.capacityOfNeighborTable_small );
  }

  fprintf(FpForLog,"\n\n");
  fflush(FpForLog);

}




void
NEIGH_displayErrorMessageForMemoryOverFlow(
										   int  iParticle
										   ,int  *numberOfNeighborParticles
										   ,int  **neighborTable
										   ,char *nameOfTable
										   ){

  int iNeigh;

  fprintf(FpForLog,"ERROR: the number of neighboring particles exceed the limit.\n");
  fprintf(FpForLog,"       Please increase the capacityOfNeighborTable.\n");
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"       %s[%d].num = %d\n",nameOfTable, iParticle, numberOfNeighborParticles[iParticle]);
  fprintf(FpForLog,"                   ");

  for(iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
    fprintf(FpForLog," %d ",neighborTable[iParticle][iNeigh] );

    if((iNeigh != 0 )&&( (iNeigh % 9) == 0)){
      fprintf(FpForLog,"\n");
      fprintf(FpForLog,"                   ");
    }

  }

  fprintf(FpForLog,"\n");


}


void
NEIGH_displayRecommendedCapacityOfNeighborTable( void ){

  int marginOfCapacity = 3;

  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "====================================================\n");
  fprintf(FpForLog, "        Recommended capacity of neighborTable       \n");
  fprintf(FpForLog, "====================================================\n");
  fprintf(FpForLog, "Recommended CapacityOfLargeNeighborTable = maximuNumberOfNeighborParticlesInLargeTable + margin\n");
  fprintf(FpForLog, "                                         = %d  + %d\n", parameter.maxArraySizeOfOccupiedNeighborTable_large, marginOfCapacity);
  fprintf(FpForLog, "                                         = %d\n", parameter.maxArraySizeOfOccupiedNeighborTable_large + marginOfCapacity);

  if( parameter.numberOfNeighborTables > 1){
    fprintf(FpForLog, "Recommended CapacityOfSmallNeighborTable = maximuNumberOfNeighborParticlesInSmallTable + margin\n");
    fprintf(FpForLog, "                                         = %d + %d\n", parameter.maxArraySizeOfOccupiedNeighborTable_small, marginOfCapacity);
    fprintf(FpForLog, "                                         = %d\n", parameter.maxArraySizeOfOccupiedNeighborTable_small + marginOfCapacity);
  }

  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "\n");
  fflush(FpForLog);

}



void
NEIGH_allocateMemoryForNeighborTable( void ){


  if(parameter.flagOfBiCG == ON){
    particle.neighborTable_large = MEMORY_allocateMemoryFor2dimIntArray( particle.totalNumber_upperLimit * 2, parameter.capacityOfNeighborTable_large, "particle.neighborTable_large" );
	particle.neighborTable_largeBoundary = MEMORY_allocateMemoryFor2dimIntArray( particle.totalNumber_upperLimit * 2, parameter.capacityOfNeighborTable_large, "particle.neighborTable_large" );
  }else{
    particle.neighborTable_large = MEMORY_allocateMemoryFor2dimIntArray( particle.totalNumber_upperLimit, parameter.capacityOfNeighborTable_large, "particle.neighborTable_large" );
	particle.neighborTable_largeBoundary = MEMORY_allocateMemoryFor2dimIntArray( particle.totalNumber_upperLimit, parameter.capacityOfNeighborTable_large, "particle.neighborTable_large" );
  }

  if( parameter.numberOfNeighborTables == 1 ) return;
  if( parameter.maxRadius_ratio == parameter.minRadius_ratio){
	parameter.numberOfNeighborTables = 1;
	return;
  }

  particle.neighborTable_small = MEMORY_allocateMemoryFor2dimIntArray( particle.totalNumber_upperLimit
                             ,parameter.capacityOfNeighborTable_small,"particle.neighborTable_small");
  particle.neighborTable_smallBoundary = MEMORY_allocateMemoryFor2dimIntArray( particle.totalNumber_upperLimit
                             ,parameter.capacityOfNeighborTable_small,"particle.neighborTable_small");

}

