#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "maxmin.h"
#include "memory.h"
#include "file.h"
#include "domain.h"
#include "other.h"
#include "bucket.h"


void
BUCKET_initializeBucket( void ){

  domain.maxNumberOfParticlesInBucket = 0;

  BUCKET_setWidthOfBucket();

  BUCKET_setNumberOfBuckets();

  BUCKET_setCapacityOfBucketAutomatically();

  BUCKET_displayInformationOfBucket();

  BUCKET_allocateMemoryForBuckets();

  BUCKET_resetParticleCounterToZero();

}



void
BUCKET_storeParticlesInBuckets( double **position ){
  
  int iParticle;
  int iX, iY, iZ;
  int flag;
  int numberOfStoredParticles;


  if(( domain.flagOfAutoSettingOfBucketCapacity == ON )&& (domain.flagOfOptimizationOfBucketMemorySize == ON )){

	BUCKET_freeMemoryOfParticleList();

	BUCKET_allocateMemoryForParticleListOfBucket( position );

  }

  BUCKET_resetParticleCounterToZero();

  for(iParticle = 0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == GHOST ) continue;

    flag = DOMAIN_checkWhetherParticleIsInDomain( iParticle, position );

    if( flag == OUT_OF_DOMAIN ){

      OTHER_changeParticleTypeIntoGhost(iParticle);
      continue;

    }

    BUCKET_findBucketWhereParticleIsStored( &iX, &iY, &iZ ,iParticle ,position );

    numberOfStoredParticles = domain.bucket[iX][iY][iZ].count;

	BUCKET_checkCapacityOfTheBucket( numberOfStoredParticles, iX, iY, iZ );

    BUCKET_storeTheParticleInBucket ( iParticle, numberOfStoredParticles, iX, iY, iZ );

  }

  BUCKET_recordMaxNumberOfParticlesInBucket();

}


void
BUCKET_countTheNumberOfParticleInEachBucket( double **position ){

  int iParticle;
  int iX, iY, iZ;
  int flag;



  for(iParticle = 0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == GHOST ) continue;

    flag = DOMAIN_checkWhetherParticleIsInDomain( iParticle, position );

    if( flag == OUT_OF_DOMAIN ) continue;

    BUCKET_findBucketWhereParticleIsStored( &iX, &iY, &iZ ,iParticle ,position );

    domain.bucket[iX][iY][iZ].count++;

  }

}




void
BUCKET_allocateMemoryForParticleListOfBucket( double **position ){

  int iX, iY, iZ;
  int capacity;


  BUCKET_resetParticleCounterToZero();

  BUCKET_countTheNumberOfParticleInEachBucket( position );

  for (iX=0 ; iX < domain.numberOfBuckets[XDIM] ; iX++) {

	for (iY=0 ; iY < domain.numberOfBuckets[YDIM] ; iY++) {

	  for (iZ=0 ; iZ < domain.numberOfBuckets[ZDIM] ; iZ++) {

		capacity = domain.bucket[iX][iY][iZ].count;

		if( capacity > 0 ){

		  domain.bucket[iX][iY][iZ].list = (int *)malloc( sizeof(int) * capacity );

		}

	  }
	}
  }

}


void
BUCKET_freeMemoryOfParticleList( void ){

  int iX, iY, iZ;

  for (iX=0 ; iX < domain.numberOfBuckets[XDIM] ; iX++) {

	for (iY=0 ; iY < domain.numberOfBuckets[YDIM] ; iY++) {

	  for (iZ=0 ; iZ < domain.numberOfBuckets[ZDIM] ; iZ++) {

		if( domain.bucket[iX][iY][iZ].count > 0 ){

		  free(domain.bucket[iX][iY][iZ].list);
		}

	  }
	}
  }
}





void
BUCKET_storeParticlesInBucketsNORMAL( double **position ){
  
  int iParticle;
  int iX, iY, iZ;
  int flag;
  int numberOfStoredParticles;

  BUCKET_resetParticleCounterToZero();


  for(iParticle = 0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == GHOST ) continue;

    flag = DOMAIN_checkWhetherParticleIsInDomain( iParticle, position );

    if( flag == OUT_OF_DOMAIN ){
      OTHER_changeParticleTypeIntoGhost(iParticle);
      continue;
    }

    BUCKET_findBucketWhereParticleIsStored( &iX, &iY, &iZ ,iParticle ,position );

    numberOfStoredParticles = domain.bucket[iX][iY][iZ].count;

    BUCKET_checkCapacityOfTheBucket( numberOfStoredParticles, iX, iY, iZ );

    BUCKET_storeTheParticleInBucket ( iParticle, numberOfStoredParticles, iX, iY, iZ );

  }


  BUCKET_recordMaxNumberOfParticlesInBucket();

}




void
BUCKET_setWidthOfBucket( void ){

  int iDim;

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    domain.bucketWidth[iDim]       = parameter.maxRadius;
    domain.bucketWidth_ratio[iDim] = parameter.maxRadius_ratio;
  }
  
}


void
BUCKET_setNumberOfBuckets( void ){

  int iDim;

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    domain.numberOfBuckets[iDim] = (int)floor(((fabs)(domain.upperLimit[iDim] - domain.lowerLimit[iDim]))/( domain.bucketWidth[iDim])) + 1;
  }

  if(NumberOfDimensions == 2){
    domain.numberOfBuckets[ZDIM] = 1;
  }

}




void
BUCKET_setCapacityOfBucketAutomatically( void ){

  int    iDim;
  double estimatedCapacity;

  if( domain.flagOfAutoSettingOfBucketCapacity == OFF) return;


  estimatedCapacity = 1.0;

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    estimatedCapacity *= (domain.bucketWidth_ratio[iDim] + 1.0);
  }

  domain.capacityOfBucket = (int)( estimatedCapacity * domain.marginRatioForSettingBucketCapacity);

  
}





void
BUCKET_displayInformationOfBucket( void ){

  int iDim;

  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "====================================================\n");
  fprintf(FpForLog, "        Buckets                                     \n");
  fprintf(FpForLog, "====================================================\n");

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    fprintf(FpForLog,"domain.bucketWidth[%s] = %lf * %lf [m] = %lf [m]\n"
	    ,FILE_returnDim(iDim), parameter.maxRadius/particle.averageDistance
            ,particle.averageDistance, domain.bucketWidth[iDim]);
  }


  for(iDim = 0; iDim < NumberOfDimensions; iDim++){
    fprintf(FpForLog,"domain.numberOfBuckets[%s] = %d\n", FILE_returnDim(iDim),domain.numberOfBuckets[iDim]);
  }


  fprintf(FpForLog,"capacityOfBucket = %d\n", domain.capacityOfBucket);

  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "\n");
  fflush(FpForLog);

}



void
BUCKET_allocateMemoryForBuckets( void ){

  int iX,iY,iZ;

  char errorMessage[256];

  domain.bucket = (structBucket***)malloc( sizeof(structBucket**) * domain.numberOfBuckets[XDIM] );

  if (domain.bucket==NULL) {
	sprintf(errorMessage,"domain.bucket is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n");
	OTHER_endProgram(errorMessage);
  }

  for (iX=0 ; iX < domain.numberOfBuckets[XDIM] ; iX++) {
	domain.bucket[iX] = (structBucket**)malloc( sizeof(structBucket*) * domain.numberOfBuckets[YDIM]);

    if (domain.bucket[iX]==NULL){
	  sprintf(errorMessage,"domain.bucket[%s] is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n"
			  ,FILE_returnDim(iX));

	  OTHER_endProgram(errorMessage);
    }

    for (iY=0 ; iY < domain.numberOfBuckets[YDIM] ; iY++) {
      domain.bucket[iX][iY] = (structBucket *)malloc( sizeof(structBucket) * domain.numberOfBuckets[ZDIM] );
      if (domain.bucket[iX][iY]==NULL){
		sprintf(errorMessage,"domain.bucket[%s][%s] is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n", FILE_returnDim(iX),FILE_returnDim(iY));

		OTHER_endProgram(errorMessage);
      }


	  if(( domain.flagOfAutoSettingOfBucketCapacity == ON )&& (domain.flagOfOptimizationOfBucketMemorySize == ON )) continue;

	  for (iZ=0 ; iZ < domain.numberOfBuckets[ZDIM] ; iZ++) {
		domain.bucket[iX][iY][iZ].list = (int *)malloc( sizeof(int) * domain.capacityOfBucket);

		if (domain.bucket[iX][iY][iZ].list == NULL){
		  sprintf(errorMessage,"domain.bucket[%s][%s][%s] is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n"
                              ,FILE_returnDim(iX), FILE_returnDim(iY), FILE_returnDim(iZ));
		  OTHER_endProgram(errorMessage);
		}
	  }

	}

  }

}


/*
void
BUCKET_allocateMemoryForBuckets( void ){

  int iX,iY,iZ;
  char errorMessage[256];

  domain.bucket = (structBucket***)malloc( sizeof(structBucket**) * domain.numberOfBuckets[XDIM] );

  if (domain.bucket==NULL) {
	sprintf(errorMessage,"domain.bucket is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n");
	OTHER_endProgram(errorMessage);
  }

  for (iX=0 ; iX < domain.numberOfBuckets[XDIM] ; iX++) {
	domain.bucket[iX] = (structBucket**)malloc( sizeof(structBucket*) * domain.numberOfBuckets[YDIM]);

    if (domain.bucket[iX]==NULL){
	  sprintf(errorMessage,"domain.bucket[%s] is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n"
			  ,FILE_returnDim(iX));

	  OTHER_endProgram(errorMessage);
    }

    for (iY=0 ; iY < domain.numberOfBuckets[YDIM] ; iY++) {
      domain.bucket[iX][iY] = (structBucket *)malloc( sizeof(structBucket) * domain.numberOfBuckets[ZDIM] );
      if (domain.bucket[iX][iY]==NULL){
		sprintf(errorMessage,"domain.bucket[%s][%s] is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n", FILE_returnDim(iX),FILE_returnDim(iY));

		OTHER_endProgram(errorMessage);
      }

	  for (iZ=0 ; iZ < domain.numberOfBuckets[ZDIM] ; iZ++) {
		domain.bucket[iX][iY][iZ].list = (int *)malloc( sizeof(int) * domain.capacityOfBucket);

		if (domain.bucket[iX][iY][iZ].list == NULL){
		  sprintf(errorMessage,"domain.bucket[%s][%s][%s] is NULL.  [in BUCKET_allocateMemoryForBuckets()]\n"
                              ,FILE_returnDim(iX), FILE_returnDim(iY), FILE_returnDim(iZ));
		  OTHER_endProgram(errorMessage);
		}
	  }

	}
  }

}
*/



void
BUCKET_checkCapacityOfTheBucket( int currentNumberOfStoredParticles, int iX, int iY, int iZ ){


  if(( domain.flagOfAutoSettingOfBucketCapacity == ON )&&( domain.flagOfOptimizationOfBucketMemorySize == ON )) return;

  if( currentNumberOfStoredParticles >= domain.capacityOfBucket){
    BUCKET_displayErrorMessageForMemoryOverFlow( iX, iY, iZ );
    OTHER_endProgram("in BUCKET_checkCapacityOfTheBucket()");
  }

}


void
BUCKET_displayErrorMessageForMemoryOverFlow( int iX, int iY, int iZ ){

  fprintf(FpForLog,"Error: The number of particles in bucket[%d][%d][%d] exceeded the limit.\n",iX, iY, iZ);
  fprintf(FpForLog,"       Please increase the capacityOfBucket.\n");
  fprintf(FpForLog,"       \n");
  fprintf(FpForLog,"       domain.bucket[%d][%d][%d].count = %d \n", iX, iY, iZ, domain.bucket[iX][iY][iZ].count);

  fflush(FpForLog);

}



void
BUCKET_storeTheParticleInBucket(int iParticle, int endOfList, int iX, int iY, int iZ){

    domain.bucket[iX][iY][iZ].list[endOfList] = iParticle;
    domain.bucket[iX][iY][iZ].count++;

}



void
BUCKET_resetParticleCounterToZero( void ){
  
  int iX, iY, iZ;  

  for(iX = 0; iX < domain.numberOfBuckets[XDIM]; iX++)
    for(iY = 0; iY < domain.numberOfBuckets[YDIM]; iY++)
      for(iZ = 0; iZ < domain.numberOfBuckets[ZDIM]; iZ++){
		domain.bucket[iX][iY][iZ].count = 0;
	  }

}



void 
BUCKET_findBucketWhereParticleIsStored( 
					   int *iX, int *iY, int *iZ
					  ,int            iParticle
					  ,double        **position
					  ){

  int place[3];
  int iDim;

  for(iDim = 0; iDim < NumberOfDimensions; iDim++){	  

   place[iDim] = (int)floor( (position[iDim][iParticle] - domain.lowerLimit[iDim])/domain.bucketWidth[iDim] );

	/*
   place[iDim] = (int)floor(((fabs)( position[iDim][iParticle] - domain.lowerLimit[iDim]))/( domain.bucketWidth[iDim]));
	*/
  }

  if( NumberOfDimensions == 2){
    place[ZDIM] = 0;
  }

  (*iX) = place[XDIM];
  (*iY) = place[YDIM];
  (*iZ) = place[ZDIM];

  BUCKET_checkErrorOfStorePlace( (*iX), (*iY), (*iZ), iParticle, position );

}




void
BUCKET_checkErrorOfStorePlace(int iX, int iY, int iZ, int iParticle, double **position){

  char errorMessage[256];
  int  errorFlag = OFF;

  if( iX <  0                            ) errorFlag = ON;
  if( iX >= domain.numberOfBuckets[XDIM]) errorFlag = ON;
  if( iY <  0                            ) errorFlag = ON;
  if( iY >= domain.numberOfBuckets[YDIM]) errorFlag = ON;

  if( NumberOfDimensions == 3){
    if( iZ <  0                            ) errorFlag = ON;
    if( iZ >= domain.numberOfBuckets[ZDIM]) errorFlag = ON;
  }

  if(errorFlag == ON){

    if(NumberOfDimensions == 3){
      sprintf(errorMessage
	    ,"Particle was not able to be stored in bucket\n [in check_error_of_storePlace()]\n iX=%d, iY=%d, iZ=%d \n position[XDIM][%d]=%lf, position[YDIM][%d]=%lf, position[ZDIM][%d]=%lf",iX,iY,iZ,iParticle, position[XDIM][iParticle],iParticle,position[YDIM][iParticle],iParticle,position[ZDIM][iParticle]);
    }else{
      sprintf(errorMessage
	    ,"Particle was not able to be stored in bucket\n [in check_error_of_storePlace()]\n iX=%d, iY=%d \n position[XDIM][%d]=%lf, position[YDIM][%d]=%lf",iX,iY,iParticle, position[XDIM][iParticle],iParticle,position[YDIM][iParticle]);
    }

    OTHER_endProgram(errorMessage);
  }

}
   


void
BUCKET_freeBuckets( void ){

  int iX, iY;

  for (iX=0 ; iX < domain.numberOfBuckets[XDIM] ; iX++) {
    for (iY=0 ; iY < domain.numberOfBuckets[YDIM] ; iY++) {
      free(domain.bucket[iX][iY]);
    }

    free(domain.bucket[iX]);
  }

  free(domain.bucket);

}





void
BUCKET_recordMaxNumberOfParticlesInBucket( void ){

  int iX, iY, iZ;  

  for(iX = 0; iX < domain.numberOfBuckets[XDIM]; iX++)
    for(iY = 0; iY < domain.numberOfBuckets[YDIM]; iY++)
      for(iZ = 0; iZ < domain.numberOfBuckets[ZDIM]; iZ++){

		if( domain.maxNumberOfParticlesInBucket < domain.bucket[iX][iY][iZ].count){
		  domain.maxNumberOfParticlesInBucket = domain.bucket[iX][iY][iZ].count;
		}

	  }
}




void
BUCKET_displayRecommendedCapacityOfBucket( void ){

  int marginOfCapacity = 3;

  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "====================================================\n");
  fprintf(FpForLog, "        Recommended capacity of Bucket              \n");
  fprintf(FpForLog, "====================================================\n");
  fprintf(FpForLog, "Recommended CapacityOfBucket = maximuNumberOfParticlesInBucket + margin\n");
  fprintf(FpForLog, "                             = %d  + %d\n", domain.maxNumberOfParticlesInBucket,  marginOfCapacity );
  fprintf(FpForLog, "                             = %d\n",       domain.maxNumberOfParticlesInBucket + marginOfCapacity );
  fprintf(FpForLog, "\n");
  fprintf(FpForLog, "\n");
  fflush(FpForLog);

}


