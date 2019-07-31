#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "mathe.h"


void
QUATERNION_multiplyQuaternionToQuaternion(   double *answer
											 ,double *q1  /*--q1: quaternion1 --*/
											 ,double *q2  /*--q2: quaternion2 --*/
											 ){
  double answer_tmp[4];
  double absoluteValue;
	

  answer_tmp[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  answer_tmp[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
  answer_tmp[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  answer_tmp[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];

  absoluteValue = (answer_tmp[0] * answer_tmp[0]) + (answer_tmp[1] * answer_tmp[1])  + (answer_tmp[2] * answer_tmp[2]) + (answer_tmp[3] * answer_tmp[3] );

  answer[0] = answer_tmp[0]/absoluteValue;
  answer[1] = answer_tmp[1]/absoluteValue;
  answer[2] = answer_tmp[2]/absoluteValue;
  answer[3] = answer_tmp[3]/absoluteValue;


}



void
QUATERNION_setRotationalMatrixUsingQuaternion(double rotationalMatrix[][3], double *quaternion)
{
  double x2 = quaternion[1] * quaternion[1];
  double y2 = quaternion[2] * quaternion[2];
  double z2 = quaternion[3] * quaternion[3];
  double xy = quaternion[1] * quaternion[2];
  double yz = quaternion[2] * quaternion[3];
  double zx = quaternion[3] * quaternion[1];
  double xw = quaternion[1] * quaternion[0];
  double yw = quaternion[2] * quaternion[0];
  double zw = quaternion[3] * quaternion[0];

  rotationalMatrix[0][0] = 1.0 - 2.0*( y2 + z2);  rotationalMatrix[0][1] = 2.0*(xy - zw);          rotationalMatrix[0][2] = 2.0 * (zx + yw);
  rotationalMatrix[1][0] = 2.0 * (xy + zw);       rotationalMatrix[1][1] = 1.0 - 2.0*(z2 + x2);  rotationalMatrix[1][2] = 2.0 * (yz - xw);
  rotationalMatrix[2][0] = 2.0 * (zx - yw);       rotationalMatrix[2][1] = 2.0 * (yz + xw);        rotationalMatrix[2][2] = 1.0 - 2.0 * (x2 + y2);

  /*
  rotationalMatrix[0][0] = 1.0 - y2 - z2;  rotationalMatrix[0][1] = xy + zw;        rotationalMatrix[0][2] = zx - yw;
  rotationalMatrix[1][0] = xy - zw;        rotationalMatrix[1][1] = 1.0 - z2 - x2;  rotationalMatrix[1][2] = yz + xw;
  rotationalMatrix[2][0] = zx + yw;        rotationalMatrix[2][1] = yz - xw;        rotationalMatrix[2][2] = 1.0 - x2 - y2;
  */

}





void
QUATERNION_resetQuaternion( double *quaternion ){

  quaternion[0] = 1.0;
  quaternion[1] = 0.0;
  quaternion[2] = 0.0;
  quaternion[3] = 0.0;


}





void
QUATERNION_makeQuaternion( double *quaternion, double *normalVector, double angle ){

  MATHE_normalizeVector( normalVector, 3);

  quaternion[0] = cos(angle/2.0);
  quaternion[1] = normalVector[XDIM] * sin(angle/2.0);
  quaternion[2] = normalVector[YDIM] * sin(angle/2.0);
  quaternion[3] = normalVector[ZDIM] * sin(angle/2.0);


}


