#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "memory.h"
#include "file.h"
#include "other.h"
#include "mathe.h"
#include "object.h"
#include "quaternion.h"
#include "forcedMotion.h"



void
FORCEDMOTION_initializeForcedMotion( structForcedMotion  *forcedMotion
				     ,char               *fileNameOfSamplingData
				     ,int                 typeNumberOfParticleToBeMoved
				     ){
  /*
  if(parameter.flagOfForcedMotionOfRigidBody != ON ) return;
  */

  FORCEDMOTION_readFileOfForcedMotion(            forcedMotion, fileNameOfSamplingData);

  FORCEDMOTION_displaySamplingDataOfForcedMotion( forcedMotion, fileNameOfSamplingData);

  FORCEDMOTION_transformUnitOfSamplingData(       forcedMotion );

  FORCEDMOTION_setInitialInterpolatedPosition(    forcedMotion );

  OBJECT_initializeObject( &forcedMotion->object, typeNumberOfParticleToBeMoved );
  /*
  OBJECT_initializeObject( &forcedMotion->object, parameter.typeNumberOfRigidParticle_forForcedMotion );
  */
}

void
FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( structForcedMotion  *forcedMotion){
    
  FILE *fp1 = NULL;
  double power;
  int  flag;
  
  power=0.0;
  flag=timer.iTimeStep-timer.iTimeStep_copy;
  
  if(timer.iTimeStep_copy!=-1){
     if(parameter.flagOfOutputOfTorqueFile == ON){
        if(flag==0){
           fp1 = FORCEDMOTION_setTorqueCalculationFile();
           power+=FORCEDMOTION_TorqueCalculationForStirredTank(1, fp1);
        }else{
           fp1 = FORCEDMOTION_setTorqueCalculationFile();
           power+=FORCEDMOTION_TorqueCalculationForStirredTank(3, fp1);
        }
     }
  }


  FORCEDMOTION_setObjectPositionAndVelocityUsingSamplingData( forcedMotion );

  OBJECT_reshapeObject(        &forcedMotion->object );
  OBJECT_setVelocityOfObject(  &forcedMotion->object );
  OBJECT_updateObjectProperty( &forcedMotion->object);
  
  if(timer.iTimeStep_copy!=-1){
     if(parameter.flagOfOutputOfTorqueFile == ON){
        if(flag==0){
           power+=FORCEDMOTION_TorqueCalculationForStirredTank(2, fp1);
           FORCEDMOTION_endTorqueCalculationFile(fp1);
           timer.iTimeStep_copy++;
        }else{
           power+=FORCEDMOTION_TorqueCalculationForStirredTank(4, fp1);
           FORCEDMOTION_endTorqueCalculationFile(fp1);
        }
     }
  }
}


void
FORCEDMOTION_setObjectPositionAndVelocityUsingSamplingData( structForcedMotion  *forcedMotion ){

  FORCEDMOTION_interpolateSamplingDataOfObjectPosition(  forcedMotion );

  FORCEDMOTION_substituteTheInterpolatedPositionForTheObjectPosition( forcedMotion );

  FORCEDMOTION_setRotationalMatrixOfObject_usingQuaternion( forcedMotion );

  /*
	FORCEDMOTION_setRotationalMatrixOfObject_usingEularAngle( forcedMotion );
  */

  FORCEDMOTION_setVelocityOfObject( forcedMotion );

}



void
FORCEDMOTION_setRotationalMatrixOfObject_usingQuaternion( structForcedMotion  *forcedMotion ){


  double angle= 0;
  double directionOfRotationalAxis[3];
  double deltaQuaternion[4];

  /*
  if(NumberOfDimensions == 3){
  */


	/*---- rotation around X axis ----*/
	angle = forcedMotion->object.rotationalPosition[XDIM];
	directionOfRotationalAxis[XDIM] = 1.0;
	directionOfRotationalAxis[YDIM] = 0.0;
	directionOfRotationalAxis[ZDIM] = 0.0;

	QUATERNION_resetQuaternion( forcedMotion->object.quaternion );

	QUATERNION_makeQuaternion( deltaQuaternion, directionOfRotationalAxis, angle );

	QUATERNION_multiplyQuaternionToQuaternion(  forcedMotion->object.quaternion
												,forcedMotion->object.quaternion
												,deltaQuaternion );


	/*---- rotation around Y axis ----*/
	angle = forcedMotion->object.rotationalPosition[YDIM];
	directionOfRotationalAxis[XDIM] = 0.0;
	directionOfRotationalAxis[YDIM] = 1.0;
	directionOfRotationalAxis[ZDIM] = 0.0;
	
	
	QUATERNION_makeQuaternion( deltaQuaternion, directionOfRotationalAxis, angle );

	QUATERNION_multiplyQuaternionToQuaternion(  forcedMotion->object.quaternion
												,forcedMotion->object.quaternion
												,deltaQuaternion );
	/*
  }
	*/



  /*---- rotation around Z axis ----*/
  angle = forcedMotion->object.rotationalPosition[ZDIM];
  directionOfRotationalAxis[XDIM] = 0.0;
  directionOfRotationalAxis[YDIM] = 0.0;
  directionOfRotationalAxis[ZDIM] = 1.0;


  QUATERNION_makeQuaternion( deltaQuaternion, directionOfRotationalAxis, angle );

  QUATERNION_multiplyQuaternionToQuaternion(  forcedMotion->object.quaternion
											 ,forcedMotion->object.quaternion
											 ,deltaQuaternion );


  QUATERNION_setRotationalMatrixUsingQuaternion( forcedMotion->object.rotationalMatrix, forcedMotion->object.quaternion);


}




void
FORCEDMOTION_setRotationalMatrixOfObject_usingEularAngle( structForcedMotion *forcedMotion ){

  double matrix_tmp[3][3];

  OBJECT_resetMatrix( forcedMotion->object.rotationalMatrix );


  if(NumberOfDimensions == 2){

	FORCEDMOTION_setMatrixOfRotationAround_zAxis( matrix_tmp, forcedMotion->object.rotationalPosition[ZDIM]);
	MATHE_multiplyMatrixByMatrix_matrixSizeIsThree(  forcedMotion->object.rotationalMatrix
													,forcedMotion->object.rotationalMatrix
													,matrix_tmp
													);


  }else if(NumberOfDimensions == 3){

	FORCEDMOTION_setMatrixOfRotationAround_xAxis( matrix_tmp, forcedMotion->object.rotationalPosition[XDIM]);
	MATHE_multiplyMatrixByMatrix_matrixSizeIsThree(  forcedMotion->object.rotationalMatrix
													,forcedMotion->object.rotationalMatrix
													,matrix_tmp
													);


	FORCEDMOTION_setMatrixOfRotationAround_yAxis( matrix_tmp, forcedMotion->object.rotationalPosition[YDIM]);
	MATHE_multiplyMatrixByMatrix_matrixSizeIsThree(  forcedMotion->object.rotationalMatrix
													,forcedMotion->object.rotationalMatrix
													,matrix_tmp
													);


	FORCEDMOTION_setMatrixOfRotationAround_zAxis( matrix_tmp, forcedMotion->object.rotationalPosition[ZDIM]);
	MATHE_multiplyMatrixByMatrix_matrixSizeIsThree(  forcedMotion->object.rotationalMatrix
													,forcedMotion->object.rotationalMatrix
													,matrix_tmp
													);
  }


}




void
FORCEDMOTION_setMatrixOfRotationAround_xAxis( double matrix[3][3], double theta ){

  matrix[0][0] = 1.0;  matrix[0][1] = 0.0;         matrix[0][2] = 0.0;
  matrix[1][0] = 0.0;  matrix[1][1] = cos(theta);  matrix[1][2] = -sin(theta);
  matrix[2][0] = 0.0;  matrix[2][1] = sin(theta);  matrix[2][2] =  cos(theta);

}



void
FORCEDMOTION_setMatrixOfRotationAround_yAxis( double matrix[3][3], double theta ){

  matrix[0][0] = cos(theta);  matrix[0][1] = 0.0;  matrix[0][2] = sin(theta);
  matrix[1][0] = 0.0;         matrix[1][1] = 1.0;  matrix[1][2] = 0.0;
  matrix[2][0] =-sin(theta);  matrix[2][1] = 0.0;  matrix[2][2] = cos(theta);

}



void
FORCEDMOTION_setMatrixOfRotationAround_zAxis( double matrix[3][3], double theta ){

  matrix[0][0] = cos(theta);  matrix[0][1] = -sin(theta);  matrix[0][2] = 0.0;
  matrix[1][0] = sin(theta);  matrix[1][1] =  cos(theta);  matrix[1][2] = 0.0;
  matrix[2][0] = 0.0;         matrix[2][1] =  0.0;         matrix[2][2] = 1.0;

}





void
FORCEDMOTION_readFileOfForcedMotion( structForcedMotion *forcedMotion
									 ,char              *fileNameOfSamplingData ){

  int  iSamplingData;
  int  iDim;
  char buf[256];
  FILE *fp;

  fp = FILE_openFile( fileNameOfSamplingData, "r" );

  fscanf(fp,"%s", buf);
  fscanf(fp,"%s", buf);

  FORCEDMOTION_scanHowToSetCenterOfRotation( fp, &forcedMotion->object.howToSetCenterOfRotation
											  ,"forcedMotion->object.howToSetCenterOfRotation");

  for(iDim=0; iDim < 3; iDim++){
	fscanf(fp,"%s", buf);
	FILE_scanDouble(  fp, &forcedMotion->object.centerOfRotation_directInput[iDim], "forcedMotion->object.centerOfRotation_directInput[iDim]");
  }

  fscanf(fp,"%s", buf);
  fscanf(fp,"%s", buf);
  FILE_scanInt(  fp, &forcedMotion->numberOfSamplingData,  "forcedMotion->numberOfSamplingData");

  fscanf(fp,"%s", buf);
  FILE_scanDouble(  fp, &forcedMotion->startingTimeInSamplingData,  "forcedMotion->startingTimeInSamplingData");

  FORCEDMOTION_allocateMemoryForStructureOfForcedMotion( forcedMotion);

  fscanf(fp,"%s", buf);

  for(iSamplingData = 0; iSamplingData < forcedMotion->numberOfSamplingData; iSamplingData++){

    FILE_scanDouble( fp, &forcedMotion->samplingTime[iSamplingData]
                       , "forcedMotion->samplingTime[iSamplingData]");


	for(iDim=0; iDim < 3; iDim++){
      FILE_scanDouble( fp, &forcedMotion->translationalPosition[iDim][iSamplingData]
                         , "forcedMotion->translationalPosition[iDim][iSamplingData]");
    }


	for(iDim=0; iDim < 3; iDim++){
	  FILE_scanDouble( fp, &forcedMotion->rotationalPosition[iDim][iSamplingData]
                          ,"forcedMotion->rotationalPosition[iDim][iSamplingData]");
	}

  }

  FILE_closeFile( fp, fileNameOfSamplingData );

}



void
FORCEDMOTION_scanHowToSetCenterOfRotation(FILE *fp, int *flag, char *variable_name){

  char buf[256];
  char errorMessage[256];

  fscanf(fp,"%s", buf);   
  if( (((strcmp( buf, "centerOfGravity" ) == 0 )
	||(strcmp( buf, "centerofgravity" ) == 0 )))
        ||((strcmp( buf, "CenterOfGravity" ) == 0 )
        ||(strcmp( buf, "CENTER_OF_GRAVITY" ) == 0 ))){

    (*flag) = CENTER_OF_GRAVITY;

  }else if( (((strcmp( buf, "centerOfWidth" ) == 0 )
	||(strcmp( buf, "centerofwidth" ) == 0 )))
        ||((strcmp( buf, "CenterOfWidth" ) == 0 )
        ||(strcmp( buf, "CENTER_OF_WIDTH" ) == 0 ))){

    (*flag) = CENTER_OF_WIDTH;

  }else if( (((strcmp( buf, "directInput" ) == 0 )
	||(strcmp( buf, "directinput" ) == 0 )))
        ||((strcmp( buf, "DirectInput" ) == 0 )
        ||(strcmp( buf, "DIRECT_INPUT" ) == 0 ))){

    (*flag) = DIRECT_INPUT;

  }else{
    sprintf(errorMessage,"Error: Parameter of \"%s\" is not adequate. [ in data-file]\n", variable_name);
    OTHER_endProgram(errorMessage);
  }

}



char*
FORCEDMOTION_returnHowToSetCenterOfRotation(int flag){

  if(flag == CENTER_OF_GRAVITY){

    return("centerOfGravity");

  }else if(flag == CENTER_OF_WIDTH){

    return("centerOfWidth");

  }else if(flag == DIRECT_INPUT){

    return("directInput");

  }else{

    return("***  ERROR  ****");

  }

}




void
FORCEDMOTION_displaySamplingDataOfForcedMotion( structForcedMotion *forcedMotion
												,char              *fileNameOfSamplingData
												){

  int  iSamplingData;
  int  iDim;


  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        %s                  \n",fileNameOfSamplingData);
  fprintf(FpForLog,"====================================================\n");
  
  fprintf(FpForLog,"#--------CenterOfRotation-----------------------------------------------------\n");
  fprintf(FpForLog,"centerOfRotation(centerOfGravity/centerOfWidth/directInput)    %s\n"
		  ,FORCEDMOTION_returnHowToSetCenterOfRotation(forcedMotion->object.howToSetCenterOfRotation));


  for(iDim=0; iDim < 3; iDim++){
	fprintf(FpForLog,"---directInput---centerOfRotation[%s](m)                     %lf\n"
			,FILE_returnDim(iDim),  forcedMotion->object.centerOfRotation_directInput[iDim]);
  }


  fprintf(FpForLog,"#--------SamplingData-------------------------------------------------------------\n");

  fprintf(FpForLog,"numberOfSamplingData       %d\n", forcedMotion->numberOfSamplingData);
  fprintf(FpForLog,"startingTime               %lf\n",forcedMotion->startingTimeInSamplingData);

  fprintf(FpForLog,"SamplingTime-----X(m)-----Y(m)-----Z(m)-----thetaX(degree)---thetaY(degree)---thetaZ(degree)\n");


  for(iSamplingData = 0; iSamplingData < forcedMotion->numberOfSamplingData; iSamplingData++){

    fprintf(FpForLog,"%lf", forcedMotion->samplingTime[iSamplingData]);

    for(iDim=0; iDim < 3; iDim++){
      fprintf(FpForLog,"\t%lf", forcedMotion->translationalPosition[iDim][iSamplingData]);
    }

    for(iDim=0; iDim < 3; iDim++){
	  fprintf(FpForLog,"\t%lf", forcedMotion->rotationalPosition[iDim][iSamplingData]);
	}


    fprintf(FpForLog,"\n");

  }

  fflush(FpForLog);

}


void
FORCEDMOTION_allocateMemoryForStructureOfForcedMotion( structForcedMotion *forcedMotion ){

  forcedMotion->samplingTime          = MEMORY_allocateMemoryFor1dimDoubleArray(forcedMotion->numberOfSamplingData,"forcedMotion->samplingTime");

  forcedMotion->translationalPosition = MEMORY_allocateMemoryFor2dimDoubleArray( 3, forcedMotion->numberOfSamplingData ,"forcedMotion->translationalPosition");

  forcedMotion->rotationalPosition    = MEMORY_allocateMemoryFor2dimDoubleArray(  3, forcedMotion->numberOfSamplingData, "forcedMotion->rotationalPosition");

}




void
FORCEDMOTION_transformUnitOfSamplingData( structForcedMotion  *forcedMotion ){

  int iDim;
  int iSamplingData;    

  for(iSamplingData = 0; iSamplingData < forcedMotion->numberOfSamplingData; iSamplingData++){

    for(iDim=0; iDim < 3; iDim++){
      
      forcedMotion->rotationalPosition[iDim][iSamplingData]  *=  ( M_PI / 180.0 ); 
      
    }

    /*
    if(NumberOfDimensions == 2){

      forcedMotion->rotationalPosition[ZDIM][iSamplingData]  *=  ( M_PI / 180.0 ); 

    }else if(NumberOfDimensions == 3){

      for(iDim=0; iDim < NumberOfDimensions; iDim++){

		forcedMotion->rotationalPosition[iDim][iSamplingData]  *=  ( M_PI / 180.0 ); 

      }

    }
    */

  }
}



void
FORCEDMOTION_setInitialInterpolatedPosition( structForcedMotion  *forcedMotion ){

  int iDim;

  FORCEDMOTION_interpolateSamplingDataOfObjectPosition(  forcedMotion );

  for(iDim=0; iDim < 3; iDim++){
	forcedMotion->initialPosition_translational[iDim] = forcedMotion->interpolatedPosition_translational[iDim];
	forcedMotion->initialPosition_rotational[iDim]    = forcedMotion->interpolatedPosition_rotational[iDim];
  }

}


void
FORCEDMOTION_interpolateSamplingDataOfObjectPosition( structForcedMotion  *forcedMotion ){

  int endOfArray;
  int iSamplingData;


  endOfArray = forcedMotion->numberOfSamplingData - 1;
  forcedMotion->correctedSimulationTime= timer.simulationTime + forcedMotion->startingTimeInSamplingData;
  /*
  if( timer.simulationTime >= forcedMotion->samplingTime[endOfArray]){
  */

  if( forcedMotion->correctedSimulationTime >= forcedMotion->samplingTime[endOfArray]){

    FORCEDMOTION_giveFinalSamplingData( forcedMotion );

  }else{

    iSamplingData = FORCEDMOTION_returnSampleNumber( forcedMotion, forcedMotion->correctedSimulationTime );

    FORCEDMOTION_interpolateTranslationalPosition( forcedMotion, iSamplingData );
    FORCEDMOTION_interpolateRotationalPosition(    forcedMotion, iSamplingData );

  }


}



void
FORCEDMOTION_giveFinalSamplingData( structForcedMotion *forcedMotion ){

  int iDim;
  int endOfArray;

  endOfArray = forcedMotion->numberOfSamplingData - 1;

  for(iDim=0; iDim < 3; iDim++){
    forcedMotion->interpolatedPosition_translational[iDim] = forcedMotion->translationalPosition[iDim][endOfArray];
  }


  for(iDim=0; iDim < 3; iDim++){
    forcedMotion->interpolatedPosition_rotational[iDim] = forcedMotion->rotationalPosition[iDim][endOfArray];
  }


  /*
  for(iDim=0; iDim < NumberOfDimensions; iDim++){
    forcedMotion->interpolatedPosition_translational[iDim] = forcedMotion->translationalPosition[iDim][endOfArray];
  }



  if( NumberOfDimensions == 2){

    forcedMotion->interpolatedPosition_rotational[ZDIM] = forcedMotion->rotationalPosition[ZDIM][endOfArray];

  }else if( NumberOfDimensions == 3){

    for(iDim=0; iDim < NumberOfDimensions; iDim++){
      forcedMotion->interpolatedPosition_rotational[iDim] = forcedMotion->rotationalPosition[iDim][endOfArray];
    }

  }else{

    OTHER_endProgram("NumberOfDimensions is not adequate.\n [in FORCEDMOTION_setObjectPositionInThePlaceOfFinalSamplingData]");

  }
  */

    
}




void
FORCEDMOTION_interpolateTranslationalPosition( structForcedMotion *forcedMotion, int iSamplingData){

  int iDim;
  /*
  for(iDim=0; iDim < NumberOfDimensions; iDim++){
  */
  for(iDim=0; iDim < 3; iDim++){

    forcedMotion->interpolatedPosition_translational[iDim] 
	  = MATHE_linearInterpolation(  forcedMotion->translationalPosition[iDim][iSamplingData    ]
								   ,forcedMotion->translationalPosition[iDim][iSamplingData - 1]
								   ,forcedMotion->correctedSimulationTime
								   ,forcedMotion->samplingTime[iSamplingData -1]
								   ,forcedMotion->samplingTime[iSamplingData]
								   ,forcedMotion->samplingTime[iSamplingData -1]
								   );

  }
}



void
FORCEDMOTION_interpolateRotationalPosition( structForcedMotion *forcedMotion, int iSamplingData ){

  int iDim;

  for(iDim=0; iDim < 3; iDim++){

    forcedMotion->interpolatedPosition_rotational[iDim] 
      = MATHE_linearInterpolation(  forcedMotion->rotationalPosition[iDim][iSamplingData    ]
				    ,forcedMotion->rotationalPosition[iDim][iSamplingData - 1]
				    ,forcedMotion->correctedSimulationTime
				    ,forcedMotion->samplingTime[iSamplingData -1]
				    ,forcedMotion->samplingTime[iSamplingData]
				    ,forcedMotion->samplingTime[iSamplingData -1]
				    );
  }

  /*
  if( NumberOfDimensions == 2 ){

    forcedMotion->interpolatedPosition_rotational[ZDIM] 
	  = MATHE_linearInterpolation(  forcedMotion->rotationalPosition[ZDIM][iSamplingData    ]
								   ,forcedMotion->rotationalPosition[ZDIM][iSamplingData - 1]
								   ,forcedMotion->correctedSimulationTime
								   ,forcedMotion->samplingTime[iSamplingData -1]
								   ,forcedMotion->samplingTime[iSamplingData]
								   ,forcedMotion->samplingTime[iSamplingData -1]
								   );
  }else{
     
    for(iDim=0; iDim < NumberOfDimensions; iDim++){

      forcedMotion->interpolatedPosition_rotational[iDim] 
		= MATHE_linearInterpolation(  forcedMotion->rotationalPosition[iDim][iSamplingData    ]
									 ,forcedMotion->rotationalPosition[iDim][iSamplingData - 1]
									 ,forcedMotion->correctedSimulationTime
									 ,forcedMotion->samplingTime[iSamplingData -1]
									 ,forcedMotion->samplingTime[iSamplingData]
									 ,forcedMotion->samplingTime[iSamplingData -1]
									 );
    }
  }
  */

}



void
FORCEDMOTION_substituteTheInterpolatedPositionForTheObjectPosition( structForcedMotion *forcedMotion ){

  int iDim;
  /*
  for(iDim=0; iDim < NumberOfDimensions; iDim++){
    forcedMotion->object.centerOfRotation[iDim] = forcedMotion->object.centerOfRotation_initial[iDim] 
	                               + (forcedMotion->interpolatedPosition_translational[iDim] 
								      - forcedMotion->initialPosition_translational[iDim] );
  }
  */

  for(iDim=0; iDim < 3; iDim++){
    forcedMotion->object.centerOfRotation[iDim] = forcedMotion->object.centerOfRotation_initial[iDim] 
	                               + (forcedMotion->interpolatedPosition_translational[iDim] 
								      - forcedMotion->initialPosition_translational[iDim] );
  }


  for(iDim=0; iDim < 3; iDim++){
    forcedMotion->object.rotationalPosition[iDim] = forcedMotion->interpolatedPosition_rotational[iDim]
		                               - forcedMotion->initialPosition_rotational[iDim];
  }

  if(NumberOfDimensions == 2){
    forcedMotion->object.rotationalPosition[XDIM] = 0.0;
    forcedMotion->object.rotationalPosition[YDIM] = 0.0;
  }

  /*
  if(NumberOfDimensions == 2){

    forcedMotion->object.rotationalPosition[XDIM] = 0.0;

    forcedMotion->object.rotationalPosition[YDIM] = 0.0;

    forcedMotion->object.rotationalPosition[ZDIM] = forcedMotion->interpolatedPosition_rotational[ZDIM] 
                                     - forcedMotion->initialPosition_rotational[ZDIM];

  }else if(NumberOfDimensions == 3){

    for(iDim=0; iDim < 3; iDim++){
      forcedMotion->object.rotationalPosition[iDim] = forcedMotion->interpolatedPosition_rotational[iDim]
		                               - forcedMotion->initialPosition_rotational[iDim];
    }

  }
  */

}




int
FORCEDMOTION_returnSampleNumber( structForcedMotion *forcedMotion, double simulationTime ){

  int iSamplingData = 0;

  for( iSamplingData = 0; iSamplingData < forcedMotion->numberOfSamplingData; iSamplingData++){
    if (forcedMotion->samplingTime[iSamplingData] > simulationTime) break;
  }

  return iSamplingData;

}


void
FORCEDMOTION_setVelocityOfObject( structForcedMotion  *forcedMotion ){

  int iDim;
  int endOfArray;
  int iSamplingData;

  double translationalDisplacement;
  double rotationalDisplacement;
  double timeDifference;

  endOfArray = forcedMotion->numberOfSamplingData - 1;
  forcedMotion->correctedSimulationTime= timer.simulationTime + forcedMotion->startingTimeInSamplingData;


  if( forcedMotion->correctedSimulationTime >= forcedMotion->samplingTime[endOfArray]){

	for( iDim=0; iDim < 3; iDim++){
	  forcedMotion->object.translationalVelocity[iDim] = 0.0;
	  forcedMotion->object.angularVelocity[iDim]       = 0.0;
	}

  }else{

    iSamplingData  = FORCEDMOTION_returnSampleNumber( forcedMotion, forcedMotion->correctedSimulationTime );
	timeDifference = ( forcedMotion->samplingTime[iSamplingData] - forcedMotion->samplingTime[iSamplingData - 1] );

	if( timeDifference != 0.0 ){

	  for( iDim=0; iDim < 3; iDim++){
		translationalDisplacement = ( forcedMotion->translationalPosition[iDim][iSamplingData ] - forcedMotion->translationalPosition[iDim][iSamplingData - 1 ] );
		rotationalDisplacement    = ( forcedMotion->rotationalPosition[   iDim][iSamplingData ] - forcedMotion->rotationalPosition[   iDim][iSamplingData - 1 ] );

		forcedMotion->object.translationalVelocity[iDim] = translationalDisplacement / timeDifference ;
		forcedMotion->object.angularVelocity[iDim]       = rotationalDisplacement    / timeDifference ;
	  }

	}else{

	  fprintf(FpForLog,"WARNING: Sampling time of forcedMotion motion is not adequate. The time difference is zero.\n");
	  fprintf(FpForLog,"         forcedMotion->samplingTime[%d] = %lf\n", iSamplingData-1, forcedMotion->samplingTime[iSamplingData-1]);
	  fprintf(FpForLog,"         forcedMotion->samplingTime[%d] = %lf\n", iSamplingData  , forcedMotion->samplingTime[iSamplingData  ]);


	  for( iDim=0; iDim < 3; iDim++){
		forcedMotion->object.translationalVelocity[iDim] = 0.0;
		forcedMotion->object.angularVelocity[iDim]       = 0.0;
	  }

	}


	if( NumberOfDimensions == 2){

	  forcedMotion->object.translationalVelocity[ZDIM] = 0.0;

	  forcedMotion->object.angularVelocity[XDIM]    = 0.0;
	  forcedMotion->object.angularVelocity[YDIM]    = 0.0;

	}

  }

}

FILE*
FORCEDMOTION_setTorqueCalculationFile(void){
    
	FILE *fp;

	// NEW!!!!! ->
	// char fileName[256];
    //
    // FORCEDMOTION_nameOfTorqueCalculationFile(fileName);
    //
	// fp=FILE_openFile(fileName,"a");

    if( timer.iTimeStep == 0 && timer.iTimeStep_copy == 0 )
      fp = FILE_openFile(parameter.nameOfOutputTorqueFile, "w");
    else
      fp = FILE_openFile(parameter.nameOfOutputTorqueFile, "a");
	// NEW!!!!! <-
    
	return(fp);
}


void
FORCEDMOTION_nameOfTorqueCalculationFile(char *fileName){
    
	char *nameOfTorqueCalculationFile;
	
	nameOfTorqueCalculationFile = "Torque.prof";
    
	strcpy(fileName, nameOfTorqueCalculationFile);
    
}

void
FORCEDMOTION_endTorqueCalculationFile(FILE *fp){
    
	char fileName[256];
    
    FORCEDMOTION_nameOfTorqueCalculationFile(fileName);
    
	FILE_closeFile(fp,fileName);
}

double
FORCEDMOTION_TorqueCalculationForStirredTank(int numberOfStatus, FILE *fp){
    
    int iParticle;
	int iDim;
    double Force[3];
    double Advection[3];
    double Torque[3];
    double Torque_previous[3];
    
    double CenterOfRotation[3];
    double PowerOfStirring;
    
    PowerOfStirring=0.0;
    
    Torque[XDIM]=0.0;
    Torque[YDIM]=0.0;
    Torque[ZDIM]=0.0;
    
#pragma omp parallel for private(iDim)
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        if(particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion){
            if(numberOfStatus == 1){
                for(iDim=0;iDim< NumberOfDimensions;iDim++){
                    particle.position1[iDim][iParticle]=particle.position[iDim][iParticle];
                    particle.velocity1[iDim][iParticle]=particle.velocity[iDim][iParticle];
                }
            }else if(numberOfStatus == 2){
                for(iDim=0;iDim< NumberOfDimensions;iDim++){
                    particle.position2[iDim][iParticle]=particle.position[iDim][iParticle];
                    particle.velocity2[iDim][iParticle]=particle.velocity[iDim][iParticle];
                }
            }else if(numberOfStatus == 3){
                for(iDim=0;iDim< NumberOfDimensions;iDim++){
                    particle.position3[iDim][iParticle]=particle.position[iDim][iParticle];
                    particle.velocity3[iDim][iParticle]=particle.velocity[iDim][iParticle];
                }
            }else if(numberOfStatus == 4){
                for(iDim=0;iDim< NumberOfDimensions;iDim++){
                    particle.position4[iDim][iParticle]=particle.position[iDim][iParticle];
                    particle.velocity4[iDim][iParticle]=particle.velocity[iDim][iParticle];
                }
            }else{
                continue;
            }
        }
    }
    
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        if(particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion){
            if(numberOfStatus == 1){
                
                CenterOfRotation[XDIM]=0.0;
                CenterOfRotation[YDIM]=particle.position_previous[YDIM][iParticle];
                CenterOfRotation[ZDIM]=0.0;
                
                Torque_previous[XDIM]=Torque[XDIM];
                Torque_previous[YDIM]=Torque[YDIM];
                Torque_previous[ZDIM]=Torque[ZDIM];
                
                for(iDim=0;iDim< NumberOfDimensions;iDim++){
                    Force[iDim]=physicalProperty.massDensity[particle.type[iParticle]]*(4.0/3.0)*M_PI*pow((particle.averageDistance/2.0),3.0)*(1.0/timer.dt)*(particle.velocity1[iDim][iParticle]-particle.velocity_previous[iDim][iParticle]);
                }
                
                Advection[XDIM]=particle.position_previous[XDIM][iParticle]-CenterOfRotation[XDIM];
                Advection[YDIM]=particle.position_previous[YDIM][iParticle]-CenterOfRotation[YDIM];
                Advection[ZDIM]=particle.position_previous[ZDIM][iParticle]-CenterOfRotation[ZDIM];
                
                MATHE_outerProduct(Torque,Advection,Force);
                
                MATHE_addVectorToVector(NumberOfDimensions,Torque,Torque,Torque_previous,1.0);
                
            }else if(numberOfStatus == 3){
                
                CenterOfRotation[XDIM]=0.0;
                CenterOfRotation[YDIM]=particle.position2[YDIM][iParticle];
                CenterOfRotation[ZDIM]=0.0;
                
                Torque_previous[XDIM]=Torque[XDIM];
                Torque_previous[YDIM]=Torque[YDIM];
                Torque_previous[ZDIM]=Torque[ZDIM];
                
                for(iDim=0;iDim< NumberOfDimensions;iDim++){
                    Force[iDim]=physicalProperty.massDensity[particle.type[iParticle]]*(4.0/3.0)*M_PI*pow((particle.averageDistance/2.0),3.0)*(1.0/timer.dt)*(particle.velocity3[iDim][iParticle]-particle.velocity2[iDim][iParticle]);
                }
                
                Advection[XDIM]=particle.position2[XDIM][iParticle]-CenterOfRotation[XDIM];
                Advection[YDIM]=particle.position2[YDIM][iParticle]-CenterOfRotation[YDIM];
                Advection[ZDIM]=particle.position2[ZDIM][iParticle]-CenterOfRotation[ZDIM];
                
                MATHE_outerProduct(Torque,Advection,Force);
                
                MATHE_addVectorToVector(NumberOfDimensions, Torque,Torque,Torque_previous,1.0);
                
            }else{
                continue;
            }
        }
    }
    PowerOfStirring+=2.0*M_PI*parameter.numberOfRotations*Torque[YDIM];
    // printf("%lf\n",PowerOfStirring);
    if (numberOfStatus <= 3)
    	fprintf(fp, "%lf ", PowerOfStirring);
    else
    	fprintf(fp, "%lf\n", PowerOfStirring);

    return(PowerOfStirring);
}



