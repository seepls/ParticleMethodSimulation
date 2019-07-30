#ifndef __FORCEDMOTION_H__
#define __FORCEDMOTION_H__


/*
void
FORCEDMOTION_initializeForcedMotion( structForcedMotion  *forcedMotion
									 ,char               *fileNameOfSamplingData
									 );
*/
void
FORCEDMOTION_initializeForcedMotion( structForcedMotion  *forcedMotion
				     ,char               *fileNameOfSamplingData
				     ,int                 typeNumberOfParticleToBeMoved
				     );

void
FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( structForcedMotion  *forcedMotion);

									 

void
FORCEDMOTION_setObjectPositionAndVelocityUsingSamplingData( structForcedMotion  *forcedMotion );



void
FORCEDMOTION_setRotationalMatrixOfObject_usingQuaternion( structForcedMotion  *forcedMotion );




void
FORCEDMOTION_setRotationalMatrixOfObject_usingEularAngle( structForcedMotion *forcedMotion );


void
FORCEDMOTION_setMatrixOfRotationAround_xAxis( double matrix[3][3], double theta );



void
FORCEDMOTION_setMatrixOfRotationAround_yAxis( double matrix[3][3], double theta );


void
FORCEDMOTION_setMatrixOfRotationAround_zAxis( double matrix[3][3], double theta );


void
FORCEDMOTION_readFileOfForcedMotion( structForcedMotion *forcedMotion
									 ,char              *fileNameOfSamplingData );



void
FORCEDMOTION_scanHowToSetCenterOfRotation(FILE *fp, int *flag, char *variable_name);



char*
FORCEDMOTION_returnHowToSetCenterOfRotation(int flag);



void
FORCEDMOTION_displaySamplingDataOfForcedMotion( structForcedMotion *forcedMotion
												,char              *fileNameOfSamplingData
												);


void
FORCEDMOTION_allocateMemoryForStructureOfForcedMotion( structForcedMotion *forcedMotion );




void
FORCEDMOTION_transformUnitOfSamplingData( structForcedMotion  *forcedMotion );



void
FORCEDMOTION_setInitialInterpolatedPosition( structForcedMotion  *forcedMotion );


void
FORCEDMOTION_interpolateSamplingDataOfObjectPosition( structForcedMotion  *forcedMotion );



void
FORCEDMOTION_giveFinalSamplingData( structForcedMotion *forcedMotion );



void
FORCEDMOTION_interpolateTranslationalPosition( structForcedMotion *forcedMotion, int iSamplingData);



void
FORCEDMOTION_interpolateRotationalPosition( structForcedMotion *forcedMotion, int iSamplingData );



void
FORCEDMOTION_substituteTheInterpolatedPositionForTheObjectPosition( structForcedMotion *forcedMotion );


int
FORCEDMOTION_returnSampleNumber( structForcedMotion *forcedMotion, double simulationTime );


void
FORCEDMOTION_setVelocityOfObject( structForcedMotion  *forcedMotion );

FILE*
FORCEDMOTION_setTorqueCalculationFile(void);

void
FORCEDMOTION_nameOfTorqueCalculationFile(char *fileName);

void
FORCEDMOTION_endTorqueCalculationFile(FILE *fp);

double
FORCEDMOTION_TorqueCalculationForStirredTank(int numberOfStatus, FILE* fp);


#endif
