//===================================================
//Program to create input files for Poisseuille Flow
//
//===================================================


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <iomanip>

#include "struct.h"
#include "file.h"
#include "init.h"
#include "geometry.h"

#define PI 3.14159265

using namespace std;


void
GEOMETRY_initialize(){
	// Set the parameters for the chosen type of geometry

	parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types;
	parameter.domainUpperLimit[2] = 0.0;
	parameter.domainLowerLimit[2] = 0.0;
	geometry.minCor[2] = 0.0;
	geometry.maxCor[2] = 0.0;
	strcpy_s (parameter.flagOfPeriodic,"off");
	strcpy_s (parameter.flagOfSingleVortex,"off");
	strcpy_s (parameter.flagOfDomainAutoSetting,"off");
	strcpy_s (parameter.flagOfInterfaceCorrection,"off");
	strcpy_s (parameter.flagOfForcedMotion,"off");
	parameter.forcedMotionParticleType = 5;
	parameter.gravity = 0.0;

	// Water Column
	if (geometry.geometryType==0){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types + 2;
		parameter.wallParticleType = parameter.number_of_Particle_Types - 2;
		parameter.dummyWallParticleType = parameter.number_of_Particle_Types - 1;
		parameter.width = 0.4;
		parameter.length = 1.0;
		geometry.dimensions = 2;
		parameter.gravity = -9.81;
		geometry.minCor[0] = 0 - 4 * parameter.l_0;
		geometry.maxCor[0] = ((int)(parameter.length / parameter.l_0))*parameter.l_0 + 4 * parameter.l_0;
		geometry.minCor[1] = 0 - 4 * parameter.l_0;
		geometry.maxCor[1] = ((int)(parameter.width / parameter.l_0))*parameter.l_0 + 2 * parameter.l_0;
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");
		parameter.domainUpperLimit[0] = parameter.length + parameter.l_0;
		parameter.domainUpperLimit[1] = 0.29;
		parameter.domainLowerLimit[0] = 0.0;
		parameter.domainLowerLimit[1] = -0.28;
	}
	// 2D Duct
	else if (((geometry.geometryType>=1) && (geometry.geometryType<=5))||(geometry.geometryType == 20)){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types + 2;
		parameter.wallParticleType = parameter.number_of_Particle_Types - 2;
		parameter.dummyWallParticleType = parameter.number_of_Particle_Types - 1;
		parameter.width = 0.5;
		parameter.length = 2.0;
		geometry.dimensions = 2;
		geometry.minCor[0] = 0;
		geometry.maxCor[0] = ((int)(parameter.length / parameter.l_0))*parameter.l_0;
		geometry.minCor[1] = -(((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0 + 4 * parameter.l_0);
		geometry.maxCor[1] = ((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0 + 4 * parameter.l_0;
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");
		parameter.domainUpperLimit[0] = parameter.length + parameter.l_0;
		parameter.domainUpperLimit[1] = 0.29;
		parameter.domainLowerLimit[0] = 0.0;
		parameter.domainLowerLimit[1] = -0.28;
		strcpy_s (parameter.flagOfPeriodic,"on");
		strcpy_s (parameter.flagOfInterfaceCorrection,"on");
		if (geometry.geometryType == 20){
			geometry.dimensions = 3;
			geometry.minCor[2] = -(((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0);
			geometry.maxCor[2] = ((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0;
		}
	}
	// Single Vortex Flow
	else if (geometry.geometryType==6){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types;
		parameter.wallParticleType = 3;
		parameter.dummyWallParticleType = 4;
		parameter.width = 2.0;
		parameter.length = 2.0;
		geometry.dimensions = 2;
		geometry.minCor[0] = 0 + parameter.l_0/2;
		geometry.maxCor[0] = 1 - parameter.l_0/2;
		geometry.minCor[1] = 0 + parameter.l_0/2;
		geometry.maxCor[1] = 1 - parameter.l_0/2;
		parameter.domainUpperLimit[0] = 1.0;
		parameter.domainUpperLimit[1] = 1.0;
		parameter.domainLowerLimit[0] = 0.0;
		parameter.domainLowerLimit[1] = 0.0;
		strcpy_s (parameter.flagOfSingleVortex,"on");
		strcpy_s (parameter.flagOfInterfaceCorrection,"on");
	}
	// 3D Bucket with Cube in Cube
	else if (geometry.geometryType==7){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types + 2;
		parameter.wallParticleType = parameter.number_of_Particle_Types - 2;
		parameter.dummyWallParticleType = parameter.number_of_Particle_Types - 1;
		geometry.dimensions = 3;
		geometry.minBucket[0] = -1;
		geometry.maxBucket[0] = 1;
		geometry.minBucket[1] = 0;
		geometry.maxBucket[1] = 1;
		geometry.minBucket[2] = -1;
		geometry.maxBucket[2] = 1;
		geometry.minCor[0] = geometry.minBucket[0] - 5 * parameter.l_0;
		geometry.maxCor[0] = geometry.maxBucket[0] + 5 * parameter.l_0;
		geometry.minCor[1] = geometry.minBucket[1] - 5 * parameter.l_0;
		geometry.maxCor[1] = geometry.maxBucket[1] + 5 * parameter.l_0;
		geometry.minCor[2] = geometry.minBucket[2] - 5 * parameter.l_0;
		geometry.maxCor[2] = geometry.maxBucket[0] + 5 * parameter.l_0;
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");

		for (int i = 0; i < 3; i++){
			geometry.center[i] = (geometry.minBucket[i] + geometry.maxBucket[i]) / 2;
			geometry.width[i] = geometry.maxBucket[i] - geometry.minBucket[i];
		}
	}
	// 3D small Cube in big Cube
	else if (geometry.geometryType==8){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types;
		parameter.wallParticleType = 3;
		parameter.dummyWallParticleType = 4;
		geometry.dimensions = 3;
		geometry.minBucket[0] = 0;
		geometry.maxBucket[0] = 1;
		geometry.minBucket[1] = 0;
		geometry.maxBucket[1] = 1;
		geometry.minBucket[2] = 0;
		geometry.maxBucket[2] = 1;
		geometry.minCor[0] = geometry.minBucket[0];
		geometry.maxCor[0] = geometry.maxBucket[0];
		geometry.minCor[1] = geometry.minBucket[1];
		geometry.maxCor[1] = geometry.maxBucket[1];
		geometry.minCor[2] = geometry.minBucket[2];
		geometry.maxCor[2] = geometry.maxBucket[0];
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");

		for (int i = 0; i < 3; i++){
			geometry.center[i] = (geometry.minBucket[i] + geometry.maxBucket[i]) / 2;
			geometry.width[i] = geometry.maxBucket[i] - geometry.minBucket[i];
			geometry.width[i] = 0.5;
		}
	}
	// 3D sphere
	else if (geometry.geometryType==9){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types;
		parameter.wallParticleType = 3;
		parameter.dummyWallParticleType = 4;
		geometry.dimensions = 3;
		geometry.minCor[0] = -0.5;
		geometry.maxCor[0] = 0.5;
		geometry.minCor[1] = -0.5;
		geometry.maxCor[1] = 0.5;
		geometry.minCor[2] = -0.5;
		geometry.maxCor[2] = 0.5;
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");
	}
	// 3D circular Poiseuille Flow
	else if (geometry.geometryType==10){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types + 2;
		parameter.wallParticleType = parameter.number_of_Particle_Types - 2;
		parameter.dummyWallParticleType = parameter.number_of_Particle_Types - 1;
		parameter.width = 0.5;
		parameter.length = 2.0;
		geometry.dimensions = 3;
		geometry.minCor[0] = 0;
		geometry.maxCor[0] = ((int)(parameter.length / parameter.l_0))*parameter.l_0;
		geometry.minCor[1] = -(((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0 + 6 * parameter.l_0);
		geometry.maxCor[1] = ((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0 + 6 * parameter.l_0;
		geometry.minCor[2] = -(((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0 + 6 * parameter.l_0);
		geometry.maxCor[2] = ((int)(parameter.width / 2 / parameter.l_0))*parameter.l_0 + 6 * parameter.l_0;
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");
		parameter.domainLowerLimit[0] = 0.0;
		parameter.domainUpperLimit[0] = parameter.length + parameter.l_0;
		parameter.domainLowerLimit[1] = -0.28;
		parameter.domainUpperLimit[1] = 0.29;
		parameter.domainLowerLimit[1] = -0.28;
		parameter.domainUpperLimit[1] = 0.29;
		strcpy_s (parameter.flagOfPeriodic,"on");
		strcpy_s (parameter.flagOfInterfaceCorrection,"on");
	}
	// 3D Stirring Tank
	else if (geometry.geometryType>=11){
		parameter.number_of_Liquid_Types = 2;
		parameter.number_of_Particle_Types = parameter.number_of_Liquid_Types + 3;
		parameter.wallParticleType = parameter.number_of_Particle_Types - 2;
		parameter.dummyWallParticleType = parameter.number_of_Particle_Types - 1;
		parameter.forcedMotionParticleType = parameter.number_of_Particle_Types -3;
		geometry.dimensions = 3;
		strcpy_s (parameter.flagOfForcedMotion,"on");
		strcpy_s (parameter.flagOfDomainAutoSetting,"on");
		geometry.diameter = 2.0;
		parameter.gravity = -9.81;
		geometry.minCor[0] = -(((int) (geometry.diameter/2/parameter.l_0))*parameter.l_0 + 6*parameter.l_0);
		geometry.maxCor[0] = ((int) (geometry.diameter/2/parameter.l_0))*parameter.l_0 + 6*parameter.l_0;
		geometry.minCor[1] = -(((int) (geometry.diameter/4/parameter.l_0))*parameter.l_0 + 6*parameter.l_0);
		geometry.maxCor[1] = 2.45;
		geometry.minCor[2] = -(((int) (geometry.diameter/2/parameter.l_0))*parameter.l_0 + 6*parameter.l_0);
		geometry.maxCor[2] = ((int) (geometry.diameter/2/parameter.l_0))*parameter.l_0 + 6*parameter.l_0;
		parameter.domainLowerLimit[0] = -3.0;
		parameter.domainUpperLimit[0] = 3;
		parameter.domainLowerLimit[1] = -7.0;
		parameter.domainUpperLimit[1] = 3;
		parameter.domainLowerLimit[1] = -3.0; // bug  ?? 
		parameter.domainUpperLimit[1] = 3; 
		strcpy_s (parameter.flagOfInterfaceCorrection,"on");
	}
}


// 0 : no particle
// 0  : fluid
// 1  : rigid
// 2  : wall (with pressure)
// 3  : wall (w/o pressure)
void GEOMETRY_particle_type(const double &x, const double &y, const double &z){
	parameter.velocity[0]=0;
	parameter.velocity[1]=0;
	parameter.velocity[2]=0;
	
	if (geometry.geometryType == 0) GEOMETRY_particle_type_column(x,y,z);

	else if (((geometry.geometryType>=1) && (geometry.geometryType<=5))||(geometry.geometryType == 20)) GEOMETRY_particle_type_duct(x,y,z);

	else if (geometry.geometryType == 6) GEOMETRY_particle_type_vortex(x,y,z);

	else if (geometry.geometryType == 7) GEOMETRY_particle_type_3D_Bucket(x,y,z);

	else if (geometry.geometryType == 8) GEOMETRY_particle_type_3D_Simple(x,y,z);

	else if (geometry.geometryType == 9) GEOMETRY_particle_type_3D_spinningBall(x,y,z);

	else if (geometry.geometryType == 10) GEOMETRY_particle_type_3D_CircularPoiseuilleFlow(x,y,z);

	else if (geometry.geometryType >= 11) GEOMETRY_particle_type_3D_StirringTank(x,y,z);

	return;
}



void GEOMETRY_particle_type_duct(const double &x, const double &y, const double &z){

    if (fabs(y) <= parameter.width / 2 + parameter.eps){
		if (x < 0 - parameter.eps || x > parameter.length + parameter.eps){
			parameter.type=-1;
			return;
		}
		else if (((geometry.geometryType == 1)||(geometry.geometryType == 20)) && (( x < 0.2 + parameter.eps) || (x >= parameter.intf - parameter.eps))){ //straight (initial)
			parameter.type=0;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - y) * (parameter.width/2 + y) / (parameter.width * parameter.width);
			return;
		}
		else if ((geometry.geometryType == 2) && (x>0.5-parameter.eps)&&(x<=1.5+parameter.eps)&&(fabs(y)<=0.2+parameter.eps)){ //rectangle
			parameter.type=0;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - y) * (parameter.width/2 + y) / (parameter.width * parameter.width);
			return;
		}
		else if ((geometry.geometryType == 3) && (x+0.5*y < 1 + parameter.eps)){ //steep diagonal
			parameter.type=0;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - y) * (parameter.width/2 + y) / (parameter.width * parameter.width);
			return;
		}
		else if ((geometry.geometryType == 4) && ((x-1)*(x-1)+y*y < 0.2*0.2)){ // circle
			parameter.type=0;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - y) * (parameter.width/2 + y) / (parameter.width * parameter.width);
			return;
		}
		else if ((geometry.geometryType == 5) && ((x-1)*(x-1)+y*y <= 0.2*0.2 + parameter.eps)&&(y>=0-parameter.eps)){ //half circle
			parameter.type=0;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - y) * (parameter.width/2 + y) / (parameter.width * parameter.width);
			return;
		}
		else {
			parameter.type=1;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - y) * (parameter.width/2 + y) / (parameter.width * parameter.width);
			return;
		}
    }
    if (fabs(y) > parameter.width / 2 - parameter.eps && fabs(y) <= parameter.width / 2 + 2 * parameter.l_0 + parameter.eps){
		parameter.type=2;
		return;
	}
    if (fabs(y) > parameter.width / 2 + 2 * parameter.l_0 - parameter.eps && fabs(y) <= parameter.width / 2 + 4 * parameter.l_0 + parameter.eps){
		parameter.type=3;
		return;
	}
	
    parameter.type=-1;
	return;
}

void GEOMETRY_particle_type_column(const double &x, const double &y, const double &z){

    if ((x >= 0)&&(x <= parameter.length + parameter.eps)&&(y >= 0)&&(y <= parameter.width)){
		if((y <= 0.5)&&(x <= 0.2)){
			parameter.type = 0;
			return;
		}
		else{
			parameter.type = -1;
			return;
		}
	}

	if((x >= 0 - 2 * parameter.l_0)&&(x <= parameter.length + 2 * parameter.l_0 + parameter.eps)&&(y >= 0 - 2 * parameter.l_0)){
		if((x >= 0)&&(x <= parameter.length + parameter.eps)&&(y >= 0)){
			parameter.type = -1;
			return;
		}
		parameter.type = parameter.wallParticleType;
		return;
	}else{
		parameter.type = parameter.dummyWallParticleType;
		return;
	}

    parameter.type=-1;
	return;
}



void GEOMETRY_particle_type_vortex(const double &x, const double &y, const double &z){

	if ((x<geometry.minCor[0]-parameter.eps)||(x>geometry.maxCor[0]+parameter.eps)||(y<geometry.minCor[1]-parameter.eps)||(y>geometry.maxCor[1]+parameter.eps)){
		parameter.type=-1;
		return;
	}

	if ((x-0.5)*(x-0.5) + (y-0.25)*(y-0.25) < 0.2*0.2 +parameter.eps){
		parameter.type=0;
		return;
	}
	else{
		parameter.type=1;
		return;
	}

}


void GEOMETRY_particle_type_3D_Bucket(const double &x, const double &y, const double &z){
	
	if ((x < geometry.minBucket[0] - 3 * parameter.l_0 - parameter.eps)||(z < geometry.minBucket[2] - 3 * parameter.l_0 - parameter.eps)||
		(x > geometry.maxBucket[0] + 3 * parameter.l_0 + parameter.eps)||(z > geometry.maxBucket[2] + 3 * parameter.l_0 + parameter.eps)||
		(y < geometry.minBucket[1] - 3 * parameter.l_0 - parameter.eps)){

			parameter.type = parameter.dummyWallParticleType;
			return;
	} else if ((x < geometry.minBucket[0] - parameter.eps)||(z < geometry.minBucket[2] - parameter.eps)||
		(x > geometry.maxBucket[0] + parameter.eps)||(z > geometry.maxBucket[2] + parameter.eps)||
		(y < geometry.minBucket[1] - parameter.eps)){

			parameter.type = parameter.wallParticleType;
			return;
	}

	if (y > geometry.maxBucket[1] + parameter.eps){
		parameter.type = -1;
		return;
	}

	if ((x > geometry.center[0] - geometry.width[0] / 4 - parameter.eps)&&(x < geometry.center[0] + geometry.width[0] / 4 + parameter.eps)&&
		(y > geometry.center[1] - geometry.width[1] / 4 - parameter.eps)&&(y < geometry.center[1] + geometry.width[1] / 4 + parameter.eps)&&
		(z > geometry.center[2] - geometry.width[2] / 4 - parameter.eps)&&(z < geometry.center[2] + geometry.width[2] / 4 + parameter.eps)){
			parameter.type = 0;
			return;
	} else{
		parameter.type = 1;
		return;
	}
}


void GEOMETRY_particle_type_3D_Simple(const double &x, const double &y, const double &z){

	if ((x > 0.25 -parameter.eps)&&(x < 0.75 + parameter.eps)&&
		(y > 0.25 -parameter.eps)&&(y < 0.75 + parameter.eps)&&
		(z > 0.25 -parameter.eps)&&(z < 0.75 + parameter.eps)){
			parameter.type = 0;
			return;
	} else{
		parameter.type = 1;
		return;
	}
}

void GEOMETRY_particle_type_3D_spinningBall(const double &x, const double &y, const double &z){
	
	double outerRadius = 0.5;
	double innerRadius = 0.25;
	
	double r = sqrt(x*x + y*y + z*z);

	if (r <= innerRadius + parameter.eps){
		parameter.type = 0;
		return;
	}else if (r <= outerRadius + parameter.eps){
		parameter.type = 1;
		return;
	}else{
		parameter.type = -1;
		return;
	}

}

void GEOMETRY_particle_type_3D_CircularPoiseuilleFlow(const double &x, const double &y, const double &z){

	double r = sqrt(y * y + z * z);

	if (fabs(r) <= parameter.width / 2 + parameter.eps){
		if (x < 0 - parameter.eps || x > parameter.length + parameter.eps){
			parameter.type=-1;
			return;
		}
		else if (( x < 0.2 + parameter.eps) || (x >= parameter.intf + parameter.eps)){ //straight (initial)
			parameter.type=0;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - r) * (parameter.width/2 + r) / (parameter.width * parameter.width);
			return;
		}
		else {
			parameter.type=1;
			parameter.velocity[0] = 4 * parameter.max_Velocity * (parameter.width/2 - r) * (parameter.width/2 + r) / (parameter.width * parameter.width);
			return;
		}
	}
	if ((fabs(r) > parameter.width / 2 - parameter.eps) && (fabs(r) <= parameter.width / 2 + 3 * parameter.l_0 + parameter.eps)){
		parameter.type=2;
		return;
	}
    if (fabs(r) > parameter.width / 2 + 2 * parameter.l_0 - parameter.eps && fabs(r) <= parameter.width / 2 + 6 * parameter.l_0 + parameter.eps){
		parameter.type=3;
		return;
	}
	
    parameter.type=-1;
	return;
}

void GEOMETRY_particle_type_3D_StirringTank(const double &x, const double &y, const double &z){

	double w = 3*parameter.l_0;  // Thickness of Wall
	double r = sqrt(x*x + z*z);
	double rd_0 = sqrt(x*x + 4*y*y + z*z);
	double rd_1 = sqrt(x*x + 4*(y+2*parameter.l_0)*(y+2*parameter.l_0) + z*z);
	double rd_2 = sqrt(x*x + 4*y*y);
	double t = 0.04; //Thickness of stirrer
	double H = 1.0; // Depth of Fluid in Tank

	if (r >= geometry.diameter/2 + 2*w - parameter.eps){
		parameter.type = -1;
		return;
	}
	// ~~~~~~~~~~Wall~~~~~~~~~~
	// -----dummyWall-----
	if ((y >= - parameter.l_0 - parameter.eps &&  r >= geometry.diameter/2 + w - parameter.eps) 
		|| (y <= -parameter.l_0 - parameter.eps && rd_1 >= geometry.diameter/2 + w - parameter.eps)){
			parameter.type = parameter.dummyWallParticleType;
			return;
	}
	// -----Wall-----
	else if ((y >= -parameter.eps && r >= geometry.diameter/2 - parameter.eps) 
		|| (y < -parameter.eps && rd_0 >= geometry.diameter/2 - parameter.eps)) {
			parameter.type = parameter.wallParticleType;
			return;
	} else{

		// ~~~~~~~~~~Stirrer~~~~~~~~~~
		// -----Paddle-----
		if (geometry.geometryType == 11){
			double d = 1.2;
			double h = 0.24;
			double hd = 0.8;

			if (r <= 2.5*t - parameter.eps && y >= -h/2){
				parameter.type = parameter.forcedMotionParticleType;
				return;
			}
			else if (abs(z) <= t + parameter.eps && abs(x) <= d/2 + parameter.eps
				&& ((y >= -h/2 + parameter.eps && y <= h/2 + parameter.eps) 
				|| (y >= hd - h/2 + parameter.eps && y <= hd + h/2 + parameter.eps))){
					parameter.type = parameter.forcedMotionParticleType;
					return;
			}
		}

		// ----- new blade -----
		else if (geometry.geometryType == 15){
			double d = 1.2;
			double h = 0.24;
			double hd = 0.8;
			double  theta ; 
			double x_position , y_position , z_position ;

			if (r <= 2.5*t - parameter.eps && y >= -h/4){
				parameter.type = parameter.forcedMotionParticleType;
				return; // rod center 
			}
			else if (abs(z) <= t + parameter.eps && abs(x) <= d/2 + parameter.eps
				&& ((y >= -h/4  + parameter.eps && y <= h/4 + parameter.eps) 
				|| (y >= hd - h/4 + parameter.eps && y <= hd + h/4 + parameter.eps) || (y >= hd/2 - h/4 + parameter.eps && y <= hd/2 +h/4 + parameter.eps) )){
					parameter.type = parameter.forcedMotionParticleType;
					return; // three arms 
			}

			for(y_position = 0.8,theta = 0  ; y_position >= 0 , theta < 2* PI ; y_position = y_position - 0.0005 , theta= theta + (0.225 * PI / 180)){ // single helical Attempt 
				//for (theta= 0 ; theta < 2* PI ; theta= theta + (0.225 * PI / 180 ) ){ // attempt for a ring at y = 0.8 : DONE 
						 x_position = d/2* cos(theta); 
						 z_position = d/2 * sin(theta);

							if( (x >= x_position - h/4 + parameter.eps && x <= x_position + h/4 + parameter.eps) &&
								(z <= z_position + h/4 + parameter.eps  && z >= z_position - h/4 + parameter.eps) && 
								(y >= y_position - h/4 + parameter.eps && y <= y_position + h/4 + parameter.eps) ){
							
									parameter.type = parameter.forcedMotionParticleType;
									return;
							}
				//}
			}


			for(y_position = 0.8,theta = 0  ; y_position >= 0 , theta < 2* PI ; y_position = y_position - 0.0005 , theta= theta + (0.225 * PI / 180)){ // double hellical attempt 
				//for (theta= 0 ; theta < 2* PI ; theta= theta + (0.225 * PI / 180 ) ){ // attempt for a ring at y = 0.8 : DONE 
						 x_position = -1*d/2* cos(theta); 
						 z_position = d/2 * sin(theta);

							if( (x >= x_position - h/4 + parameter.eps && x <= x_position + h/4 + parameter.eps) &&
								(z <= z_position + h/4 + parameter.eps  && z >= z_position - h/4 + parameter.eps) && 
								(y >= y_position - h/4 + parameter.eps && y <= y_position + h/4 + parameter.eps) ){
							
									parameter.type = parameter.forcedMotionParticleType;
									return;
							}
				//}
			}

			 /*
			for ( y_position = 0.8 ; y_position >= 0 ; y_position= y_position - 0.0005){
				for (theta= 0 ; theta < 2* PI ; theta= theta + (0.225 * PI / 180 ) ){
						 x_position = d/2* cos(theta); 
						 z_position = d/2 * sin(theta);

							if( (x >= x_position - h/4 && x <= x_position + h/4) &&
								(z <= z_position && z >= z_position) && (y >= y_position - h/4 + parameter.eps && y <= y_position + h/4 + parameter.eps)){
							
									parameter.type = parameter.forcedMotionParticleType;
									return;
							}
				}

			}
			*/
		
		}
			
		
		
		// -----Anchor-----
		else if(geometry.geometryType == 12){
			double d = 1.6;
			double h = 0.8;
			double a = 0.16;
			double del_y = a;
			double rd_3 = sqrt(x*x + 4*(y-del_y)*(y-del_y));
			double offset = a - 0.5*sqrt(a*d - a*a); // for smootheness

			if (r <= 2.5*t - parameter.eps && y >= -h/2) {
				parameter.type = parameter.forcedMotionParticleType;
				return;
			} else if (abs(z) <= t + parameter.eps) {  // toge
				// curve
				if (y <= - parameter.eps && rd_3 >= d/2 + parameter.eps && rd_2 <= d/2 + parameter.eps){
					parameter.type = parameter.forcedMotionParticleType;
					return;
				}
				// straight part on the side
				if (y >= - parameter.eps && y <= h/2 && abs(x) >= d/2 - a + parameter.eps && abs(x) <= d/2 + parameter.eps){
					parameter.type = parameter.forcedMotionParticleType;
					return;
				}
				// for smoothness
				if (y >= offset && y <= - parameter.eps && abs(x) >= d/2 - a + parameter.eps && abs(x) <= d/2 - a/2 + parameter.eps){
					parameter.type = parameter.forcedMotionParticleType;
					return;
				}
			}
		}
		// -----Maxblend-----
		else if(geometry.geometryType == 13){
			double d = 1.4;
			double h = 2.0;
			double a = 0.14;
			double b = 1.68;
			bool ret = false;
			
			if (r <= 2.5*t - parameter.eps && y >= -d/4) {  // jiku
				parameter.type = parameter.forcedMotionParticleType;
				return;
			} else if (abs(z) <= t + parameter.eps) {
				// is it inside the "AND symbol"?
				if (y >= -parameter.eps && y <= h + parameter.eps && abs(x) <= d/2 + parameter.eps)
					ret = true;
				if (y <= -parameter.eps && rd_2 <= d/2 + parameter.eps)
					ret = true;
				// hole
				if (y >= (h - b)/2 - parameter.eps && y <= (h + b)/2 - parameter.eps) {
					if (abs(x) >= 0.1 && abs(x) <= 0.1 + 4*0.04)
						ret = false;
					else if (abs(x) >= 0.1 + 7*0.04 && abs(x) <= 0.1 + 11*0.04)
						ret = false;
				}
				if (ret){
					parameter.type = parameter.forcedMotionParticleType;
					return;
				}
			}
		/*
		// ------- DFouble Helical -----
 		else if (geometry.geometryType == 14 ){
 			double r = 0.8 // helical ring distance from center 
 			double pitch = 1.0  // set up pitch length one complete round : here set to the whole height
 			double connectingRodRadius = // rods connecting the helical ribbon to the center rod . 
 			double centerRodRadius = // radius of center rod  // central rotating shaft 
 			double h
 			
 		}
 		*/


		


 	}
		// ~~~~~~~~~~Fluid~~~~~~~~~~
		if (y <= 0.3 * H && parameter.number_of_Particle_Types == 6){
			parameter.type = 2;
			return;
		}else if (y <= 0.7 * H && parameter.number_of_Particle_Types == 6){
			parameter.type = 1;
			return;
		}else if (y <= 0.3 * H && parameter.number_of_Particle_Types == 5){
			parameter.type = 1;
			return;
		}else if (y <= H){
			parameter.type = 0;
			return;
		}
	}

	// ~~~~~~~~~~Ghost~~~~~~~~~~
	parameter.type = -1;
	return;
}
