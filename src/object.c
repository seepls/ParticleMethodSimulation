#include <stdio.h>
#include <math.h>
#include "define.h"
#include "extern.h"
#include "struct.h"
#include "memory.h"
#include "file.h"
#include "copy.h"
#include "mathe.h"
#include "quaternion.h"
#include "other.h"
#include "object.h"



void
OBJECT_initializeObject( structObject *object, int particleType ){

  object->particleType = particleType;

  OBJECT_countNumberOfComponents(          object );

  OBJECT_allocateMemoryForObjectStructure( object );

  OBJECT_setNumberListOfComponents(        object );

  OBJECT_calculateCenterOfGravity(         object );

  OBJECT_calculateCenterOfWidth(           object );

  OBJECT_setCenterOfRotation(              object );

  OBJECT_memorizeInitialCenterOfRotation(  object );

  OBJECT_setRelativePosition(              object );

  QUATERNION_resetQuaternion(              object->quaternion);

  OBJECT_displayStateOfObject(             object );

  OBJECT_updateObjectProperty(             object );

}





void
OBJECT_countNumberOfComponents( structObject *object ){

  int  iParticle;
  int  flagOfComponent;
  char errorMessage[256];

  object->numberOfComponents = 0;

  for( iParticle = 0;  iParticle < particle.totalNumber;  iParticle++ ){
    
    flagOfComponent = NO;

    flagOfComponent = OBJECT_checkWhetherTheParticleIsComponent( iParticle, object );

    if( flagOfComponent == YES ){
      object->numberOfComponents++;
    }
  }

  if(object->numberOfComponents == 0){

	sprintf(errorMessage,"The number of object particles( rigidParticle or moveWallParticle)  is zero. [in OBJECT_countNumberOfComponents()]\n Please make the rigidBodyCalculation(or sloshingCalculation) off in %s.\n",parameter.nameOfDataFile);

    OTHER_endProgram( errorMessage );

  }

}




int 
OBJECT_checkWhetherTheParticleIsComponent( int iParticle , structObject *object ){


  int flagOfComponent = NO;

  if( particle.type[iParticle] == object->particleType) flagOfComponent = YES;

  if( (object->particleType == parameter.wallType)||(object->particleType == parameter.dummyWallType)){

    if( particle.type[iParticle] == parameter.wallType     ) flagOfComponent = YES;
    if( particle.type[iParticle] == parameter.dummyWallType) flagOfComponent = YES;

  }

  return flagOfComponent;

}



void
OBJECT_allocateMemoryForObjectStructure( structObject *object ){

  if( object->numberOfComponents > 0 ){
	object->numberListOfComponents        = MEMORY_allocateMemoryFor1dimIntArray(       object->numberOfComponents, "object->numberListOfComponents");
	object->relativePosition              = MEMORY_allocateMemoryFor2dimDoubleArray( 3, object->numberOfComponents, "object->relativePosition");
	object->relativePosition_rotated      = MEMORY_allocateMemoryFor2dimDoubleArray( 3, object->numberOfComponents, "object->relativePosition_rotated");
	object->surfaceFlag                   = MEMORY_allocateMemoryFor1dimIntArray(       object->numberOfComponents, "object->surfaceFlag");
  }

}



void
OBJECT_setNumberListOfComponents( structObject *object ){

  int iParticle;
  int iComponent;
  int flagOfComponent;

  iComponent = 0;

  for( iParticle = 0;  iParticle < particle.totalNumber;  iParticle++ ){

    flagOfComponent = OBJECT_checkWhetherTheParticleIsComponent( iParticle, object );

    if( flagOfComponent == YES ){
	  
	  OBJCECT_addTheParticleInListOfComponents( iParticle, object, iComponent );
	  iComponent++;

    }
  }

}


void
OBJCECT_addTheParticleInListOfComponents( int iParticle, structObject *object, int iComponent ){

  object->numberListOfComponents[ iComponent ] = iParticle; 

}




void
OBJECT_calculateCenterOfGravity( structObject *object ){

  int iParticle;
  int iComponent;

  int iDim;

  for(iDim=0; iDim < 3; iDim++){
    object->centerOfGravity[iDim] = 0.0;
  }

  for(iComponent = 0; iComponent < object->numberOfComponents;  iComponent++){
    iParticle = object->numberListOfComponents[iComponent];

    for(iDim=0; iDim < 3; iDim++){
      object->centerOfGravity[iDim] += particle.position[iDim][iParticle];
    }
  }


  for(iDim=0; iDim < 3; iDim++){
    object->centerOfGravity[iDim] /= (double)( object->numberOfComponents);
  }

  /*
  for(iDim=0; iDim < NumberOfDimensions; iDim++){
    object->centerOfGravity[iDim] = 0.0;
  }

  for(iComponent = 0; iComponent < object->numberOfComponents;  iComponent++){
    iParticle = object->numberListOfComponents[iComponent];

    for(iDim=0; iDim < NumberOfDimensions; iDim++){
      object->centerOfGravity[iDim] += particle.position[iDim][iParticle];
    }
  }


  for(iDim=0; iDim < NumberOfDimensions; iDim++){
    object->centerOfGravity[iDim] /= (double)( object->numberOfComponents);
  }
  */


}





void
OBJECT_calculateCenterOfWidth( structObject *object ){


  int iDim;
  int iParticle;
  int iComponent;


  OBJECT_setInitialMaxMinPositionOfComponents(  object->maxPositionOfComponents
											   ,object->minPositionOfComponents
											   ,object
											   );


  for(iComponent = 1; iComponent < object->numberOfComponents;  iComponent++){
    iParticle = object->numberListOfComponents[iComponent];
    /*
    for(iDim=0; iDim < NumberOfDimensions; iDim++){
    */
    for(iDim=0; iDim < 3; iDim++){

	  if( object->maxPositionOfComponents[iDim] < particle.position[iDim][iParticle]){
		object->maxPositionOfComponents[iDim] = particle.position[iDim][iParticle];
	  }

	  if( object->minPositionOfComponents[iDim] > particle.position[iDim][iParticle]){
		object->minPositionOfComponents[iDim] = particle.position[iDim][iParticle];
	  }

    }
  }

  /*
  for(iDim=0; iDim < NumberOfDimensions; iDim++){
  */
  for(iDim=0; iDim < 3; iDim++){
    object->centerOfWidth[iDim] = 0.5 * (object->maxPositionOfComponents[iDim] + object->minPositionOfComponents[iDim]);
  }


}





void
OBJECT_setInitialMaxMinPositionOfComponents(  double *maxPositionOfComponents
											 ,double *minPositionOfComponents
											 ,structObject   *object ){

  int iDim;
  int iParticle;
  int iComponent;

  iComponent = 0;  

  /*
  for(iDim=0; iDim < NumberOfDimensions; iDim++){
  */
  for(iDim=0; iDim < 3; iDim++){
    iParticle = object->numberListOfComponents[iComponent];      
	maxPositionOfComponents[iDim] = particle.position[iDim][iParticle];
	minPositionOfComponents[iDim] = particle.position[iDim][iParticle];
  }

}




void
OBJECT_setCenterOfRotation( structObject *object ){

  int iDim;
  /*
  if( object->howToSetCenterOfRotation == CENTER_OF_GRAVITY){

	for(iDim=0; iDim < NumberOfDimensions; iDim++){
	  object->centerOfRotation[iDim] = object->centerOfGravity[iDim];
	}

  }else if( object->howToSetCenterOfRotation == CENTER_OF_WIDTH){

	for(iDim=0; iDim < NumberOfDimensions; iDim++){
	  object->centerOfRotation[iDim] = object->centerOfWidth[iDim];
	}

  }else if( object->howToSetCenterOfRotation == DIRECT_INPUT){

	for(iDim=0; iDim < NumberOfDimensions; iDim++){
	  object->centerOfRotation[iDim] = object->centerOfRotation_directInput[iDim];
	}
  */
  if( object->howToSetCenterOfRotation == CENTER_OF_GRAVITY){

	for(iDim=0; iDim < 3; iDim++){
	  object->centerOfRotation[iDim] = object->centerOfGravity[iDim];
	}

  }else if( object->howToSetCenterOfRotation == CENTER_OF_WIDTH){

	for(iDim=0; iDim < 3; iDim++){
	  object->centerOfRotation[iDim] = object->centerOfWidth[iDim];
	}

  }else if( object->howToSetCenterOfRotation == DIRECT_INPUT){

	for(iDim=0; iDim < 3; iDim++){
	  object->centerOfRotation[iDim] = object->centerOfRotation_directInput[iDim];
	}

  }else{
	fprintf(stderr,  "howToSetCenterOfRotation is wrong.\n");
	fprintf(FpForLog,"howToSetCenterOfRotation is wrong.\n");
	fflush(FpForLog);
  }

}






void
OBJECT_memorizeInitialCenterOfRotation( structObject *object ){


  COPY_copy1dimDoubleArray( 3
			   ,object->centerOfRotation_initial
			   ,object->centerOfRotation          );

}





void
OBJECT_displayStateOfObject( structObject *object ){

  int iDim;

  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"       State Of Object                              \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"particleType                              = %d\n",object->particleType);
  fprintf(FpForLog,"The total number of components(particles) = %d\n",object->numberOfComponents);
  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"centerOfGravity[%s]                       = %e[m]\n",FILE_returnDim(iDim), object->centerOfGravity[iDim]);
  }
  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"centerOfWidth[%s]                         = %e[m]\n",FILE_returnDim(iDim), object->centerOfWidth[iDim]);
  }
  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"centerOfRotation[%s]                      = %e[m]\n",FILE_returnDim(iDim), object->centerOfRotation[iDim]);
  }
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");
  fflush(FpForLog);

}



void
OBJECT_setRelativePosition( structObject *object ){

  int iParticle;
  int iComponent;
  int iDim;

  for(iComponent = 0; iComponent < object->numberOfComponents; iComponent++){
    iParticle = object->numberListOfComponents[iComponent];
    /*
    for(iDim=0; iDim < NumberOfDimensions; iDim++){
    */
    for(iDim=0; iDim < 3; iDim++){
      object->relativePosition[iDim][iComponent]         = particle.position[iDim][iParticle] - object->centerOfRotation[iDim];
      object->relativePosition_rotated[iDim][iComponent] = particle.position[iDim][iParticle] - object->centerOfRotation[iDim];
    }

    if( NumberOfDimensions == 2 ){
      object->relativePosition[ZDIM][iComponent] = 0.0;
      object->relativePosition_rotated[ZDIM][iComponent] = 0.0;
    }

  }
}





void
OBJECT_reshapeObject( structObject *object ){

  int iParticle;
  int iComponent;
  int iDim;


  double rotatedRelativePosition[3];



  for(iComponent = 0; iComponent < object->numberOfComponents; iComponent++){

    OBJECT_initializeRotatedRelativePosition( rotatedRelativePosition, iComponent, object );


    OBJECT_rotateRelativePosition( rotatedRelativePosition, object );

    iParticle = object->numberListOfComponents[iComponent];
    /*
    for(iDim=0; iDim < NumberOfDimensions; iDim++){
    */
    for(iDim=0; iDim < 3; iDim++){
      particle.position[iDim][iParticle] = rotatedRelativePosition[iDim] + object->centerOfRotation[iDim];
    }

    /*
    for(iDim=0; iDim < NumberOfDimensions; iDim++){
    */
    for(iDim=0; iDim < 3; iDim++){
	  object->relativePosition_rotated[iDim][iComponent] = rotatedRelativePosition[iDim];
    }

  }


  OBJECT_writeObjectPositionInFile( object );


}



void
OBJECT_initializeRotatedRelativePosition(
										 double       *rotatedRelativePosition
										 ,int          iComponent
										 ,structObject *object
										 ){

  rotatedRelativePosition[XDIM] = object->relativePosition[XDIM][iComponent];
  rotatedRelativePosition[YDIM] = object->relativePosition[YDIM][iComponent];
  rotatedRelativePosition[ZDIM] = object->relativePosition[ZDIM][iComponent];

}




void
OBJECT_rotateRelativePosition(
							  double       *rotatedRelativePosition
							  ,structObject *object
							  ){


  MATHE_linearTransform( rotatedRelativePosition
						  ,object->rotationalMatrix
						  ,rotatedRelativePosition );



}


void
OBJECT_resetMatrix( double matrix[3][3] ){

  matrix[0][0] = 1.0;  matrix[0][1] = 0.0;  matrix[0][2] = 0.0;
  matrix[1][0] = 0.0;  matrix[1][1] = 1.0;  matrix[1][2] = 0.0;
  matrix[2][0] = 0.0;  matrix[2][1] = 0.0;  matrix[2][2] = 1.0;

}




void
OBJECT_setVelocityOfObject( structObject *object ){

  int iComponent, iParticle;

  for(iComponent = 0; iComponent < object->numberOfComponents; iComponent++){

    iParticle = object->numberListOfComponents[iComponent];


    if(NumberOfDimensions == 2){

	  particle.velocity[XDIM][iParticle] = object->translationalVelocity[XDIM] 
                                     - (object->angularVelocity[ZDIM] * object->relativePosition_rotated[YDIM][iComponent]);

	  particle.velocity[YDIM][iParticle] = object->translationalVelocity[YDIM] 
		                             + (object->angularVelocity[ZDIM] * object->relativePosition_rotated[XDIM][iComponent]);


	  particle.velocity[ZDIM][iParticle] = 0.0;


    }else if(NumberOfDimensions == 3){


	  particle.velocity[XDIM][iParticle] = object->translationalVelocity[XDIM] 
		                             + object->angularVelocity[YDIM] * object->relativePosition_rotated[ZDIM][iComponent] 
		                             - object->angularVelocity[ZDIM] * object->relativePosition_rotated[YDIM][iComponent];

	  particle.velocity[YDIM][iParticle] = object->translationalVelocity[YDIM] 
		                             + object->angularVelocity[ZDIM] * object->relativePosition_rotated[XDIM][iComponent] 
		                             - object->angularVelocity[XDIM] * object->relativePosition_rotated[ZDIM][iComponent];

	  particle.velocity[ZDIM][iParticle] = object->translationalVelocity[ZDIM] 
		                             + object->angularVelocity[XDIM] * object->relativePosition_rotated[YDIM][iComponent] 
		                             - object->angularVelocity[YDIM] * object->relativePosition_rotated[XDIM][iComponent];


	}
  }

}



void
OBJECT_updateObjectProperty( structObject *object ){

  int iDim; 

  for( iDim=0; iDim < 3; iDim++){
	object->centerOfGravity_previous[iDim] = object->centerOfGravity[iDim];
  }

  for( iDim=0; iDim < 3; iDim++){
	object->centerOfRotation_previous[iDim] = object->centerOfRotation[iDim]; 
  }


}




void
OBJECT_writeObjectPositionInFile( structObject *object ){

  static FILE *fp = NULL;

  double rotationalPosition[3];


  if( fp == NULL ){
	fp = FILE_openFile( "output.objectMotion", "w");
	fprintf( fp, "time\tcenterOfRotation[XDIM]\tcenterOfRotation[YDIM]\tcenterOfRotation[ZDIM]\tangleOfRotation[XDIM]\tangleOfRotation[YDIM]\tangleOfRotation[ZDIM]\tquaternion[0]\tquaternion[1]\tquaternion[2]\tquaternion[3]\n");
	/*
	fprintf( fp, "time\tcenterOfRotation[XDIM]\tcenterOfRotation[YDIM]\tcenterOfRotation[ZDIM]\tangleOfRotation[XDIM]\tangleOfRotation[YDIM]\tangleOfRotation[ZDIM]\n");
	*/
  }

  OBJECT_calculateRotationalPosition( rotationalPosition, object->quaternion );

  fprintf( fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
		   ,timer.simulationTime
		   ,object->centerOfRotation[XDIM]
		   ,object->centerOfRotation[YDIM]
		   ,object->centerOfRotation[ZDIM]
		   ,rotationalPosition[XDIM]
		   ,rotationalPosition[YDIM]
		   ,rotationalPosition[ZDIM]
		   ,object->quaternion[0]
		   ,object->quaternion[1]
		   ,object->quaternion[2]
		   ,object->quaternion[3]
		   );



  fflush(fp);


}



void
OBJECT_calculateRotationalPosition( double *rotationalPosition, double *quaternion){

  int    iAxis;
  double innerProduct;

  double rotationalMatrix[3][3];
  double directionVector[3][3];
  double directionVector_rotated[3][3];
  double directionVector_tmp[3][3];
  double radian;


  directionVector[0][XDIM] = 1.0;
  directionVector[0][YDIM] = 0.0;
  directionVector[0][ZDIM] = 0.0;

  directionVector[1][XDIM] = 0.0;
  directionVector[1][YDIM] = 1.0;
  directionVector[1][ZDIM] = 0.0;

  directionVector[2][XDIM] = 0.0;
  directionVector[2][YDIM] = 0.0;
  directionVector[2][ZDIM] = 1.0;


  QUATERNION_setRotationalMatrixUsingQuaternion( rotationalMatrix, quaternion );

  
  for( iAxis=0; iAxis < 3; iAxis++){
	MATHE_linearTransform( directionVector_rotated[iAxis]
						   ,rotationalMatrix
						   ,directionVector[iAxis]);
  }


  directionVector_tmp[0][XDIM] = 0.0;
  directionVector_tmp[0][YDIM] = directionVector_rotated[1][YDIM];
  directionVector_tmp[0][ZDIM] = directionVector_rotated[1][ZDIM];

  directionVector_tmp[1][XDIM] = directionVector_rotated[0][XDIM];
  directionVector_tmp[1][YDIM] = 0.0;
  directionVector_tmp[1][ZDIM] = directionVector_rotated[0][ZDIM];

  directionVector_tmp[2][XDIM] = directionVector_rotated[0][XDIM];
  directionVector_tmp[2][YDIM] = directionVector_rotated[0][YDIM];
  directionVector_tmp[2][ZDIM] = 0.0;


  for( iAxis=0; iAxis < 3; iAxis++){
	MATHE_normalizeVector( directionVector_tmp[iAxis], 3);
  }


  innerProduct =  MATHE_innerProductOfVectors( directionVector_tmp[2], directionVector[0], 3);
  radian = acos( innerProduct );
  rotationalPosition[2] = (radian / M_PI ) * 180.0;

  if( directionVector_rotated[0][YDIM] < 0.0 ){
	rotationalPosition[2] *= (-1.0);
  }



  innerProduct =  MATHE_innerProductOfVectors( directionVector_tmp[1], directionVector[0], 3);
  radian = acos( innerProduct );
  rotationalPosition[1] = (radian / M_PI ) * 180.0;

  if( directionVector_rotated[2][XDIM] < 0.0 ){
	rotationalPosition[1] *= (-1.0);
  }




  innerProduct =  MATHE_innerProductOfVectors( directionVector_tmp[0], directionVector[1], 3);
  radian = acos( innerProduct );
  rotationalPosition[0] = (radian / M_PI ) * 180.0;

  if( directionVector_rotated[2][YDIM] > 0.0 ){
	rotationalPosition[0] *= (-1.0);
  }



  /*
  innerProduct =  MATHE_innerProductOfVectors( directionVector_tmp[2], directionVector[0], 3);
  radian = acos( innerProduct );
  rotationalPosition[2] = (radian / M_PI ) * 180.0;


  innerProduct =  MATHE_innerProductOfVectors( directionVector_tmp[1], directionVector[0], 3);
  radian = acos( innerProduct );
  rotationalPosition[1] = (radian / M_PI ) * 180.0;


  innerProduct =  MATHE_innerProductOfVectors( directionVector_tmp[0], directionVector[1], 3);
  radian = acos( innerProduct );
  rotationalPosition[0] = (radian / M_PI ) * 180.0;
  */


}

