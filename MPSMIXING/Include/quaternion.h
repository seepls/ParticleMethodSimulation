void
QUATERNION_multiplyQuaternionToQuaternion(   double *answer
											 ,double *q1  /*--q1: quaternion1 --*/
											 ,double *q2  /*--q2: quaternion2 --*/
											 );

void
QUATERNION_setRotationalMatrixUsingQuaternion(double rotationalMatrix[][3], double *quaternion);


void
QUATERNION_resetQuaternion( double *quaternion );


void
QUATERNION_makeQuaternion( double *quaternion, double *normalVector, double angle );

