#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "weight.h"
#include "extern.h"
#include "interface.h"
#include "buoyancy.h"
#include "distance.h"
#include "file.h"
#include "gravity.h"
#include "neigh.h"
#include "other.h"
#include "memory.h"
#include "mathe.h"
#include "inverse.h"
#include "density.h"


void
INTERFACE_CalculateInterfaceArea(void){

	int i,j;

	if (parameter.flagOfInterfaceCalculation == OFF) return;

	if ((fmod(timer.simulationTime + TOO_SMALL_VALUE,timer.intervalTimeOfInterface)> timer.dt)&&(timer.dt != 0)) return;
	
	if (parameter.flagOfInterfaceCompare == ON){
		INTERFACE_calculateInterfaceCompare();
		return;
	}

	for(i=0;i<parameter.numberOfParticleTypes;i++){

		if (i == GHOST) continue;
		if (i == parameter.dummyWallType) continue;
		if (i == parameter.wallType) continue;

		parameter.totalInterfaceArea[i] = 0;
		parameter.totalLiquidInterfaceArea[i] = 0;

		for (j=0;j<parameter.numberOfParticleTypes;j++){

			if (j == GHOST) continue;
			if (j == parameter.dummyWallType) continue;
			if (j == parameter.wallType) continue;
			if (j == i) continue;

			INTERFACE_CalculateInterfaceAreaBetweenIAndJ(i,j);

			parameter.totalInterfaceArea[i] += parameter.interfaceArea[i][j];
			if (j == parameter.wallType) continue;
			if (j == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;
			parameter.totalLiquidInterfaceArea[i] += parameter.interfaceArea[i][j];
		}
	}


	INTERFACE_writeInterfaceFile(0,1);
}



void
INTERFACE_calculateInterfaceCompare(void){

	FILE *fp;

	int i;

	fp = FILE_openFile("Interface.area", "a");

	fprintf(fp, "\n");
	fprintf(fp, "%lf", timer.simulationTime);

	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);

	// for Single Vortex Flow
	if (parameter.flagOfSingleVortex == ON){
		parameter.flagOfInterfaceCorrection = OFF;
		parameter.flagOfFabs = OFF;
		parameter.flagOfCombinedInterfaceCalculation = OFF;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);
		INTERFACE_CalculateInterfaceAreaBetweenIAndJ(1, 0);
		INTERFACE_writeInterfaceFileCompare(1, 0);
		parameter.flagOfCombinedInterfaceCalculation = ON;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);

		parameter.flagOfInterfaceCorrection = ON;
		parameter.flagOfFabs = OFF;
		parameter.flagOfCombinedInterfaceCalculation = OFF;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);
		INTERFACE_CalculateInterfaceAreaBetweenIAndJ(1, 0);
		INTERFACE_writeInterfaceFileCompare(1, 0);
		parameter.flagOfCombinedInterfaceCalculation = ON;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);

		parameter.orderOfLsmps = 1;
		parameter.flagOfCombinedInterfaceCalculation = OFF;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);
		INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(1, 0);
		INTERFACE_writeInterfaceFileCompare(1, 0);
		parameter.flagOfCombinedInterfaceCalculation = ON;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);

		parameter.orderOfLsmps = 2;
		parameter.flagOfCombinedInterfaceCalculation = OFF;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);
		INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(1, 0);
		INTERFACE_writeInterfaceFileCompare(1, 0);
		parameter.flagOfCombinedInterfaceCalculation = ON;
		INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
		INTERFACE_writeInterfaceFileCompare(0, 1);
		return;
	}


	// for Wall and unCombined
	if(1){
	parameter.flagOfInterfaceCorrection = ON;
	parameter.orderOfLsmps = 1;
	parameter.flagOfFabs = OFF;
	parameter.flagOfCombinedInterfaceCalculation = OFF;
	/*INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(parameter.wallType,0);
	INTERFACE_writeInterfaceFileCompare(parameter.wallType,0);
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(parameter.wallType,1);
	INTERFACE_writeInterfaceFileCompare(parameter.wallType,1);*/
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(parameter.wallType,10);
	INTERFACE_writeInterfaceFileCompare(parameter.wallType,10);
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(10,parameter.wallType);
	INTERFACE_writeInterfaceFileCompare(10,parameter.wallType);

	parameter.flagOfCombinedInterfaceCalculation = ON;
	/*INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(parameter.wallType,0);
	INTERFACE_writeInterfaceFileCompare(parameter.wallType,0);
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(parameter.wallType,1);
	INTERFACE_writeInterfaceFileCompare(parameter.wallType,1);*/
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(parameter.wallType,10);
	INTERFACE_writeInterfaceFileCompare(parameter.wallType,10);
	//INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(10,parameter.wallType);
	//INTERFACE_writeInterfaceFileCompare(10,parameter.wallType);
	return;}

	//CorrMatrix
	if(0){
	parameter.flagOfInterfaceCorrection = OFF;
	parameter.flagOfFabs = OFF;
	parameter.flagOfCombinedInterfaceCalculation = ON;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0,1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	parameter.flagOfInterfaceCorrection = ON;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0,1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	return;}

	//Conventional Gradient
	parameter.flagOfCombinedInterfaceCalculation = ON;
	parameter.flagOfInterfaceCorrection = OFF;
	parameter.flagOfFabs = OFF;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	//Absolute distance
	parameter.flagOfFabs = ON;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);

	//Corrected Gradient
	parameter.flagOfInterfaceCorrection = ON;
	parameter.flagOfFabs = OFF;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	//Absolute distance
	parameter.flagOfFabs = ON;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJ(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	//return;

	//LSMPS 1st order
	parameter.orderOfLsmps = 1;
	parameter.flagOfFabs = OFF;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	//Absolute distance
	parameter.flagOfFabs = ON;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	return;

	//LSMPS 2nd order
	parameter.orderOfLsmps = 2;
	parameter.flagOfFabs = OFF;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
	//INTERFACE_CalculateInterfaceAreaBetweenIAndJSplit(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);
	//Absolute distance
	/*parameter.flagOfFabs = ON;
	INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(0, 1);
	INTERFACE_writeInterfaceFileCompare(0, 1);*/
	return;

}


void
INTERFACE_CalculateInterfaceAreaBetweenIAndJ(int iType, int jType){

	int iParticle, jParticle, iNeigh, iDim;

	double xji;
	double yji;
	double zji = 0.0;
	double d_cx;
	double d_cy;
	double d_cz = 0.0;

	double distanceIJ;
	double distanceIJ_squared;
	double ConcentrationGradient;
	double weight;

	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	double tot;
	double nablaF[MAX_DIM];

	double corrMatrixIJ[MAX_DIM][MAX_DIM];

	double interfaceArea = 0.0;

	int    localFlagOfFabs;

	parameter.interfaceArea[iType][jType]=0;

	parameter.radiusOfLaplacianForInterface = parameter.radiusOfGradient;

	NEIGH_selectNeighborTable(   &numberOfNeighborParticles
		,&neighborTable
		,parameter.radiusOfLaplacianForInterface
		,&neighborTablePeriodic
		);

	// Loop through all the particles
#pragma omp parallel for private(jParticle,iNeigh,iDim,xji,yji,zji,d_cx,d_cy,d_cz,distanceIJ,distanceIJ_squared,ConcentrationGradient,weight,tot,nablaF,corrMatrixIJ,localFlagOfFabs) reduction(+:interfaceArea)
	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

		// set iParticle to desired particle Type
		if (iType != 10){
			if (parameter.flagOfCombinedInterfaceCalculation == ON){
				if (jType != 10)
					if ((particle.type[iParticle] != iType)&&(particle.type[iParticle] != jType)) continue;}
			else
				if (particle.type[iParticle] != iType) continue;
		}
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
			
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaF[iDim] = 0.0;
		}

		// check if absolute value should be used for the calculation
		localFlagOfFabs = OFF;
		if (parameter.flagOfFabs == ON)
			localFlagOfFabs = INTERFACE_setFabs(iParticle, numberOfNeighborParticles[iParticle], neighborTable[iParticle], neighborTablePeriodic[iParticle], iType, jType);

		// calculate correction Matrix
		if (parameter.flagOfInterfaceCorrection == ON)
			INTERFACE_calculateCorrectionMatrix(iParticle, numberOfNeighborParticles[iParticle], neighborTable[iParticle], neighborTablePeriodic[iParticle], corrMatrixIJ);

		ConcentrationGradient = 0.0;

		// Loop over the neighbor particles
		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			// set jParticle to desired particle type
			if (jType != 10){
				if (parameter.flagOfCombinedInterfaceCalculation == ON){
					if ((particle.type[iParticle] != jType)&&(particle.type[jParticle] != jType)) continue;
					if (iType != 10)
						if ((particle.type[jParticle] != iType)&&(particle.type[jParticle] != jType)) continue;}
				else
					if (particle.type[jParticle] != jType) continue;
			}
			if ((iType != 10)&&(parameter.flagOfCombinedInterfaceCalculation == ON))
				if ((particle.type[iParticle] != iType)&&(particle.type[jParticle] != iType)) continue;
			if (particle.type[jParticle] == particle.type[iParticle]) continue;
			if (particle.type[jParticle] == parameter.dummyWallType) continue;

			// Calculate distance between i and j
			xji = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
			if((neighborTablePeriodic[iParticle][iNeigh]==1)&&(parameter.flagOfPeriodicX==ON)){
				xji += parameter.pipeLength;
			}
			else if((neighborTablePeriodic[iParticle][iNeigh]==2)&&(parameter.flagOfPeriodicX==ON)){
				xji -= parameter.pipeLength;
			}
			yji = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);

			distanceIJ_squared = xji*xji + yji*yji;

			if(NumberOfDimensions == 3){
				zji = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
				distanceIJ_squared += zji*zji;
			}


			distanceIJ = sqrt(distanceIJ_squared);

			// calculate weight of jParticle
			weight = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForInterface);

			// Calculate Gradient
			ConcentrationGradient = weight / distanceIJ_squared;

			if (localFlagOfFabs == ON){
				//xji = fabs(xji);
				yji = fabs(yji);
				if(NumberOfDimensions == 3){
					//zji = fabs(yji);
				}
			}

			if(parameter.flagOfInterfaceCorrection == OFF){
				d_cx = ConcentrationGradient * xji;
				d_cy = ConcentrationGradient * yji;
				if(NumberOfDimensions == 3){
					d_cz = ConcentrationGradient * zji;
				}
			} else if (parameter.flagOfInterfaceCorrection == ON){
				d_cx = ConcentrationGradient * (corrMatrixIJ[XDIM][XDIM] * xji + corrMatrixIJ[XDIM][YDIM] * yji + corrMatrixIJ[XDIM][ZDIM] * zji);
				d_cy = ConcentrationGradient * (corrMatrixIJ[YDIM][XDIM] * xji + corrMatrixIJ[YDIM][YDIM] * yji + corrMatrixIJ[YDIM][ZDIM] * zji);
				if(NumberOfDimensions == 3){
					d_cz = ConcentrationGradient * (corrMatrixIJ[ZDIM][XDIM] * xji + corrMatrixIJ[ZDIM][YDIM] * yji + corrMatrixIJ[ZDIM][ZDIM] * zji);
				}
			}

			// Add jParticles share to iParticles gradient
			nablaF[XDIM] += d_cx;
			nablaF[YDIM] += d_cy;

			if(NumberOfDimensions == 3){
				nablaF[ZDIM] += d_cz;
			}
		}

		// Multiply constants to the gradient
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaF[iDim] = nablaF[iDim] * NumberOfDimensions * particle.averageDistance * particle.averageDistance / parameter.nZeroOfGradient;
			if(NumberOfDimensions == 3)
				nablaF[iDim] = nablaF[iDim] * particle.averageDistance;
		}

		// optional: use a Threshold
		if (NumberOfDimensions == 2)
			INTERFACE_Threshold(nablaF,0,particle.averageDistance / 2, particle.averageDistance / 4);
		else if (NumberOfDimensions == 3)
			INTERFACE_Threshold(nablaF,0,particle.averageDistance * particle.averageDistance / 2, particle.averageDistance * particle.averageDistance / 4);

		if (parameter.flagOfCombinedInterfaceCalculation == OFF){
			for (iDim = 0; iDim < NumberOfDimensions; iDim++)
				nablaF[iDim] = nablaF[iDim] * 2;
		}

		// Calculate Eucledian length and add it to the total interface area
		tot=0;
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			tot += nablaF[iDim] * nablaF[iDim];
		interfaceArea += sqrt(tot);

		//Store Particles Gradient for VTK File
		//if (((iType==0)||(iType==1))&&((jType==0)||(jType==1)))
			for (iDim = 0; iDim < NumberOfDimensions; iDim++)
				particle.interfaceGradient[iDim][iParticle] = nablaF[iDim];

	}

	// Save interfaceArea between iType and jType
	parameter.interfaceArea[iType][jType] = interfaceArea;
}


void 
INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(int iType, int jType){

	int iParticle, jParticle, iNeigh, iDim, jDim;

	double xji;
	double yji;
	double zji = 0.0;

	double distanceIJ;
	double distanceIJ_squared;
	double weight;

	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	double tot;
	double nablaF[MAX_DIM];

	double Mi[10][10];
	double bi[10];
	double Hrs[10][10];

	double rs = particle.averageDistance;// 0.8 * parameter.radiusOfLaplacianForInterface;

	double interfaceArea = 0.0;

	int    localFlagOfFabs;

	int    matrixSize;

	double p[10];

	bi[2] = 0;

	// set size of LSMPS matrices and vectors
	if (NumberOfDimensions == 2){
		if (parameter.orderOfLsmps == 1)
			matrixSize = 2;
		else if (parameter.orderOfLsmps == 2)
			matrixSize = 5;
	}else if (NumberOfDimensions == 3){
		if (parameter.orderOfLsmps == 1)
			matrixSize = 3;
		else if (parameter.orderOfLsmps == 2)
			matrixSize = 9;
	}

	
	parameter.interfaceArea[iType][jType]=0;
	
	parameter.radiusOfLaplacianForInterface = parameter.radiusOfGradient;

	NEIGH_selectNeighborTable(   &numberOfNeighborParticles
		,&neighborTable
		,parameter.radiusOfLaplacianForInterface
		,&neighborTablePeriodic
		);

	// loop through all the particles
#pragma omp parallel for private(jParticle,iNeigh,iDim,jDim,xji,yji,zji,distanceIJ,distanceIJ_squared,weight,tot,nablaF,Mi,bi,Hrs,localFlagOfFabs,p) reduction(+:interfaceArea)
	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
		
		// set iParticle to desired particle Type
		if (iType != 10){
			if (parameter.flagOfCombinedInterfaceCalculation == ON){
				if (jType != 10)
					if ((particle.type[iParticle] != iType)&&(particle.type[iParticle] != jType)) continue;}
			else
				if (particle.type[iParticle] != iType) continue;
		}
		if ((iType == 10)&&(parameter.flagOfCombinedInterfaceCalculation == OFF)&&(particle.type[iParticle] == jType)) continue;
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
			
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaF[iDim] = 0.0;
		}

		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaF[iDim] = 0.0;
		}

		// check if absolute value should be used for the calculation
		localFlagOfFabs = OFF;
		if (parameter.flagOfFabs == ON)
			//localFlagOfFabs = ON;
			localFlagOfFabs = INTERFACE_setFabs(iParticle, numberOfNeighborParticles[iParticle], neighborTable[iParticle], neighborTablePeriodic[iParticle], iType, jType);

		// initialize LSMPS vectors and matrices
		for(iDim=0; iDim < matrixSize; iDim++){
			bi[iDim] = 0.0;
			for(jDim=0; jDim < matrixSize; jDim++){
				Mi[iDim][jDim] = 0.0;
				if(iDim==jDim)
					Hrs[iDim][jDim] = 1/rs;
				else
					Hrs[iDim][jDim] = 0.0;
			}
		}

		// Loop over the neighbor particles
		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			// set jParticle to desired particle type
			//if (jType != 10){
				//if (parameter.flagOfCombinedInterfaceCalculation == ON){
					//if ((particle.type[iParticle] != jType)&&(particle.type[jParticle] != jType)) continue;
					//if (iType != 10){
						//if ((particle.type[jParticle] != iType)&&(particle.type[jParticle] != jType)) continue;}
					//else if ((particle.type[iParticle] != jType)&&(particle.type[jParticle] != jType)) continue;
				//}
				//else
					//if (particle.type[jParticle] != jType) continue;
			//}
			//if ((iType != 10)&&(parameter.flagOfCombinedInterfaceCalculation == ON))
				//if ((particle.type[iParticle] != iType)&&(particle.type[jParticle] != iType)) continue;
			//if (particle.type[jParticle] == particle.type[iParticle]) continue;
			//if (particle.type[jParticle] == parameter.dummyWallType) continue;

			// calculate distances
			xji = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
			if((neighborTablePeriodic[iParticle][iNeigh]==1)&&(parameter.flagOfPeriodicX==ON)){
				xji += parameter.pipeLength;
			}
			else if((neighborTablePeriodic[iParticle][iNeigh]==2)&&(parameter.flagOfPeriodicX==ON)){
				xji -= parameter.pipeLength;
			}
			yji = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);
			distanceIJ_squared = xji*xji + yji*yji;
			if(NumberOfDimensions == 3){
				zji = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
				distanceIJ_squared += zji*zji;
			}
			distanceIJ = sqrt(distanceIJ_squared);

			// calculate weight
			weight = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForInterface);

			//setting vector p
			if (NumberOfDimensions == 2){
				p[0] = xji / rs;
				p[1] = yji / rs;
				if (parameter.orderOfLsmps == 2){
					p[2] = xji * xji / (rs * rs);
					p[3] = xji * yji / (rs * rs);
					p[4] = yji * yji / (rs * rs);
				}
			}
			else if (NumberOfDimensions == 3){
				p[0] = xji / rs;
				p[1] = yji / rs;
				p[2] = zji / rs;
				if (parameter.orderOfLsmps == 2){
					p[3] = xji * xji / (rs * rs);
					p[4] = yji * yji / (rs * rs);
					p[5] = zji * zji / (rs * rs);
					p[6] = xji * yji / (rs * rs);
					p[7] = yji * zji / (rs * rs);
					p[8] = xji * zji / (rs * rs);
				}
			}
			
			//calculating vector bi and Matrix Mi
			for (iDim = 0; iDim < matrixSize; iDim++){
				if ((particle.type[jParticle] != particle.type[iParticle])&&
					(particle.type[iParticle] != parameter.dummyWallType)&&(particle.type[jParticle] != parameter.dummyWallType)){
						if ((iType==10)||((iType != 10)&&((particle.type[iParticle] == iType)||(particle.type[jParticle] == iType)))){
							if ((jType == 10)||((jType != 10)&&((particle.type[iParticle] == jType)||(particle.type[jParticle] == jType)))){
								/*if (localFlagOfFabs == OFF)
									bi[iDim] += weight * p[iDim];
								else if (localFlagOfFabs == ON)
									bi[iDim] += weight * fabs(p[iDim]);*/
								if ((localFlagOfFabs == ON)&(iDim == YDIM))
									bi[iDim] += weight * fabs(p[iDim]);
								else
									bi[iDim] += weight * p[iDim];
							}
						}
				}
				for (jDim = 0; jDim < matrixSize; jDim++){
					Mi[iDim][jDim] += weight * p[iDim] * p[jDim];
				}
			}
		}

		// Invert Matrix Mi
		if((bi[0] != 0)||(bi[1] != 0)||(bi[2] != 0))
			INVERSE_matrixInverse(Mi,matrixSize);

		// Calculate gradient per dimension
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaF[iDim] = 0;
			for (jDim = 0; jDim < matrixSize; jDim++)
				nablaF[iDim] += Mi[iDim][jDim] * bi[jDim];
			nablaF[iDim] *= Hrs[iDim][iDim];
		}

		// multiply with constants
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaF[iDim] = nablaF[iDim] * particle.averageDistance * particle.averageDistance;
			if(NumberOfDimensions == 3)
				nablaF[iDim] = nablaF[iDim] * particle.averageDistance;
		}

		// optional: use a Threshold
		if (NumberOfDimensions == 2)
			INTERFACE_Threshold(nablaF,0,particle.averageDistance / 2, particle.averageDistance / 4);
		else if (NumberOfDimensions == 3)
			INTERFACE_Threshold(nablaF,0,particle.averageDistance * particle.averageDistance / 2, particle.averageDistance * particle.averageDistance / 4);

		if (parameter.flagOfCombinedInterfaceCalculation == OFF){
			for (iDim = 0; iDim < NumberOfDimensions; iDim++)
				nablaF[iDim] = nablaF[iDim] * 2;
		}

		//Add Euclidean Length of gradient of particle i to the total interface Area
		tot=0;
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			tot += nablaF[iDim] * nablaF[iDim];

		interfaceArea += sqrt(tot);

		//Store Particles Gradient for VTK File
		if (((iType==0)||(iType==1))&&((jType==0)||(jType==1)))
			for (iDim = 0; iDim < NumberOfDimensions; iDim++)
				particle.interfaceGradient[iDim][iParticle] = nablaF[iDim];
	}

	// Save interfaceArea between iType and jType
	parameter.interfaceArea[iType][jType] = interfaceArea;
}

void 
INTERFACE_CalculateInterfaceAreaBetweenIAndJSplit(int iType, int jType){

	int iParticle, jParticle, iNeigh, iDim, jDim;

	double xji;
	double yji;
	double zji = 0.0;
	double rji[MAX_DIM];
	double normVector[MAX_DIM];
	double normBi[MAX_DIM];
	double normMi[10][10];
	double normLength;
	double vectorProduct;
	double upWeight;
	
	double distanceIJ;
	double distanceIJ_squared;
	double weight;

	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	double particleInterfaceUp;
	double particleInterfaceDown;
	//double nablaF[MAX_DIM];
	double nablaFUp[MAX_DIM];
	double nablaFDown[MAX_DIM];

	double Mi[10][10];
	//double bi[10];
	double biUp[10];
	double biDown[10];
	double Hrs[10][10];

	double rs = particle.averageDistance;// 0.8 * parameter.radiusOfLaplacianForInterface;

	double interfaceArea = 0.0;

	int    localFlagOfFabs;

	int    matrixSize;

	double p[10];

	DENSITY_calculateParticleNumberDensity(particle.position);//remove later

	biUp[2] = 0;
	biDown[2] = 0;

	if (NumberOfDimensions == 2){
		if (parameter.orderOfLsmps == 1)
			matrixSize = 2;
		else if (parameter.orderOfLsmps == 2)
			matrixSize = 5;
	}else if (NumberOfDimensions == 3){
		if (parameter.orderOfLsmps == 1)
			matrixSize = 3;
		else if (parameter.orderOfLsmps == 2)
			matrixSize = 9;
	}


	parameter.interfaceArea[iType][jType]=0;

	parameter.radiusOfLaplacianForInterface = parameter.radiusOfGradient;

	NEIGH_selectNeighborTable(   &numberOfNeighborParticles
		,&neighborTable
		,parameter.radiusOfLaplacianForInterface
		,&neighborTablePeriodic
		);
//#pragma omp parallel for private(jParticle,iNeigh,iDim,jDim,xji,yji,zji,distanceIJ,distanceIJ_squared,weight,tot,nablaF,Mi,bi,Hrs,localFlagOfFabs,p) reduction(+:interfaceArea)
	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

		if ((particle.type[iParticle] != iType)&&(particle.type[iParticle] != jType)) continue;
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaFUp[iDim] = 0.0;
			nablaFDown[iDim] = 0.0;
			normBi[iDim] = 0.0;
		}

		localFlagOfFabs = OFF;
		if (parameter.flagOfFabs == ON)
			localFlagOfFabs = INTERFACE_setFabs(iParticle, numberOfNeighborParticles[iParticle], neighborTable[iParticle], neighborTablePeriodic[iParticle], iType, jType);

		for(iDim=0; iDim < matrixSize; iDim++){
			biUp[iDim] = 0.0;
			biDown[iDim] = 0.0;
			for(jDim=0; jDim < matrixSize; jDim++){
				Mi[iDim][jDim] = 0.0;
				normMi[iDim][jDim] = 0.0;
				if(iDim==jDim)
					Hrs[iDim][jDim] = 1/rs;
				else
					Hrs[iDim][jDim] = 0.0;
			}
		}


		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];
			if ((particle.type[jParticle] != iType)&&(particle.type[jParticle] != jType)&&(particle.type[jParticle] == particle.type[iParticle])) continue;

			//Distances
			rji[XDIM] = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
			if((neighborTablePeriodic[iParticle][iNeigh]==1)&&(parameter.flagOfPeriodicX==ON)){
				rji[XDIM] += parameter.pipeLength;
			}
			else if((neighborTablePeriodic[iParticle][iNeigh]==2)&&(parameter.flagOfPeriodicX==ON)){
				rji[XDIM] -= parameter.pipeLength;
			}
			rji[YDIM] = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);
			distanceIJ_squared = rji[XDIM]*rji[XDIM] + rji[YDIM]*rji[YDIM];
			if(NumberOfDimensions == 3){
				rji[ZDIM] = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
				distanceIJ_squared += rji[ZDIM]*rji[ZDIM];
			}
			distanceIJ = sqrt(distanceIJ_squared);

			//weight
			weight = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForInterface);

			for (iDim = 0; iDim < NumberOfDimensions; iDim++){
				normBi[iDim] += weight * rji[iDim] / (distanceIJ * particle.particleNumberDensity_eachParticleType[iParticle][jType]);
				for (jDim=0; jDim < matrixSize; jDim++)
					normMi[iDim][jDim] += weight * rji[iDim] * rji[jDim];
			}
		}

		//Invert Mi
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			for (jDim=0; jDim < matrixSize; jDim++)
				normMi[iDim][jDim] *= NumberOfDimensions / parameter.nZeroOfGradient;
		//INVERSE_matrixInverse(normMi,NumberOfDimensions);

		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			normVector[iDim] = normBi[iDim];//0;
			//for (jDim = 0; jDim < matrixSize; jDim++)
				//normVector[iDim] += normMi[iDim][jDim] * normBi[jDim];
		}

		//Normalize the normal Vector
		normLength = 0;
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			normLength += normVector[iDim] * normVector[iDim];
		normLength = sqrt(normLength);
		if(normLength > 0)
			printf("normLength: %f\n", normLength);
		/*
		if ((particle.type[iParticle] == iType)&&(particle.particleNumberDensity_eachParticleType[iParticle][jType]!=0)){
			upWeight = normLength / particle.particleNumberDensity_eachParticleType[iParticle][jType];
			printf("normLength: %f, particleNumberDensity: %f, upWeight: %f\n", normLength,particle.particleNumberDensity_eachParticleType[iParticle][jType],upWeight);
		}
		if ((particle.type[iParticle] == jType)&&(particle.particleNumberDensity_eachParticleType[iParticle][iType]!=0)){
			upWeight = normLength / particle.particleNumberDensity_eachParticleType[iParticle][iType];
			printf("normLength: %f, particleNumberDensity: %f, upWeight: %f", normLength,particle.particleNumberDensity_eachParticleType[iParticle][iType],upWeight);
		}*/
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			normVector[iDim] = normVector[iDim] / normLength;



		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			if ((particle.type[jParticle] != iType)&&(particle.type[jParticle] != jType)) continue;

			//Distances
			xji = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
			if((neighborTablePeriodic[iParticle][iNeigh]==1)&&(parameter.flagOfPeriodicX==ON)){
				xji += parameter.pipeLength;
			}
			else if((neighborTablePeriodic[iParticle][iNeigh]==2)&&(parameter.flagOfPeriodicX==ON)){
				xji -= parameter.pipeLength;
			}
			yji = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);
			distanceIJ_squared = xji*xji + yji*yji;
			if(NumberOfDimensions == 3){
				zji = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
				distanceIJ_squared += zji*zji;
			}
			distanceIJ = sqrt(distanceIJ_squared);

			//weight
			weight = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForInterface);

			//setting vector p
			if (NumberOfDimensions == 2){
				p[0] = xji / rs;
				p[1] = yji / rs;
				if (parameter.orderOfLsmps == 2){
					p[2] = xji * xji / (rs * rs);
					p[3] = xji * yji / (rs * rs);
					p[4] = yji * yji / (rs * rs);
				}
			}
			else if (NumberOfDimensions == 3){
				p[0] = xji / rs;
				p[1] = yji / rs;
				p[2] = zji / rs;
				if (parameter.orderOfLsmps == 2){
					p[3] = xji * xji / (rs * rs);
					p[4] = yji * yji / (rs * rs);
					p[5] = zji * zji / (rs * rs);
					p[6] = xji * yji / (rs * rs);
					p[7] = yji * zji / (rs * rs);
					p[8] = xji * zji / (rs * rs);
				}
			}

			//check normal Vector
			rji[XDIM] = xji;
			rji[YDIM] = yji;
			rji[ZDIM] = zji;
			vectorProduct = MATHE_innerProductOfVectors(rji,normBi,NumberOfDimensions);
			
			//calculating vector bi and Matrix Mi
			for (iDim = 0; iDim < matrixSize; iDim++){
				if (particle.type[jParticle] != particle.type[iParticle]){
					if ((vectorProduct >= 0)||(normLength < particle.averageDistance * 0.01))
						biUp[iDim] += weight * p[iDim];
					else
						biDown[iDim] += weight * p[iDim];
				}
				for (jDim = 0; jDim < matrixSize; jDim++){
					Mi[iDim][jDim] += weight * p[iDim] * p[jDim];
				}
			}
		}

		//Invert Matrix Mi
		if((biUp[0] != 0)||(biUp[1] != 0)||(biUp[2] != 0)||(biDown[0] != 0)||(biDown[1] != 0)||(biDown[2] != 0))
			INVERSE_matrixInverse(Mi,matrixSize);

		//Calculate gradient per dimension
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaFUp[iDim] = 0;
			nablaFDown[iDim] = 0;
			for (jDim = 0; jDim < matrixSize; jDim++){
				nablaFUp[iDim] += Mi[iDim][jDim] * biUp[jDim];
				nablaFDown[iDim] += Mi[iDim][jDim] * biDown[jDim];
			}
			nablaFUp[iDim] *= Hrs[iDim][iDim];
			nablaFDown[iDim] *= Hrs[iDim][iDim];
		}

		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			nablaFUp[iDim] = nablaFUp[iDim] * particle.averageDistance * particle.averageDistance;
			nablaFDown[iDim] = nablaFDown[iDim] * particle.averageDistance * particle.averageDistance;
			if(NumberOfDimensions == 3){
				nablaFUp[iDim] = nablaFUp[iDim] * particle.averageDistance;
				nablaFDown[iDim] = nablaFDown[iDim] * particle.averageDistance;
			}
		}

		//Threshold
		//INTERFACE_Threshold(nablaF,0,particle.averageDistance / 2, particle.averageDistance / 4);

		//Add Euclidean Length of gradient of particle i to the total interface Area
		particleInterfaceUp=0;
		particleInterfaceDown=0;
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			particleInterfaceUp += nablaFUp[iDim] * nablaFUp[iDim];
			particleInterfaceDown += nablaFDown[iDim] * nablaFDown[iDim];
		}

		interfaceArea += sqrt(particleInterfaceUp);
		interfaceArea += sqrt(particleInterfaceDown);

		//Store Particles Gradient for VTK File
		if (((iType==0)||(iType==1))&&((jType==0)||(jType==1)))
			for (iDim = 0; iDim < NumberOfDimensions; iDim++)
				particle.interfaceGradient[iDim][iParticle] = fabs(nablaFUp[iDim]) + fabs(nablaFDown[iDim]);
	}

	
	parameter.interfaceArea[iType][jType] = interfaceArea;
}

void
INTERFACE_Threshold(double nablaF[MAX_DIM], double lowerThreshold, double upperThreshold, double criteria){

	int iDim;
	double length;
	
	length = 0;

	for (iDim = 0; iDim < NumberOfDimensions; iDim++)
		length += nablaF[iDim] * nablaF[iDim];
	length = sqrt(length);

	if(length < criteria)
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			nablaF[iDim] = 0;
	else if (length >= criteria){
		for (iDim = 0; iDim < NumberOfDimensions; iDim++)
			nablaF[iDim] = nablaF[iDim] * upperThreshold / length;
		//length = upperThreshold;
	}

	

}


int
INTERFACE_setFabs(int iParticle, int numberOfNeighborParticles, int *neighborTable, int *neighborTablePeriodic, int iType, int jType){
	int iNeigh, jParticle;
	double distanceIJ, weight;
	double totalWeight = 0;
	double otherWeight = 0;
	int result;

	result = OFF;

	for (iNeigh=0; iNeigh <  numberOfNeighborParticles; iNeigh++){
		jParticle = neighborTable[iNeigh];
		//if ((particle.type[jParticle] == iType)||(particle.type[jParticle] == jType)){
			distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iNeigh]);
			weight = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForInterface);
			totalWeight += weight;
			if ((particle.type[iParticle]!=particle.type[jParticle])&&((particle.type[jParticle] == iType)||(particle.type[jParticle] == jType)))
				otherWeight += weight;
		//}
	}
	if((otherWeight > 0.55 * totalWeight)&&(NumberOfDimensions == 2))
		result = ON;
	else if((NumberOfDimensions == 3)&&(otherWeight > 0.5 * totalWeight))
		result = ON;
	return result;
}



void
INTERFACE_calculateCorrectionMatrix(int iParticle, int numberOfNeighborParticles, int *neighborTable, int *neighborTablePeriodic, double corrMatrixIJ[MAX_DIM][MAX_DIM]){

	int iNeigh, i, j, jParticle;
	double xji,yji,zji;
	double statVolume;

	double distanceIJ;
    double distanceIJ_squared;
    double weight;

	double matrix[MAX_DIM][MAX_DIM];

	statVolume=0;
	for(i=0; i < MAX_DIM; i++){
		for(j=0; j < MAX_DIM; j++){
			corrMatrixIJ[i][j] = 0.0;
			matrix[i][j] = 0.0;
		}
	}

	for(iNeigh=0; iNeigh <  numberOfNeighborParticles; iNeigh++){
		jParticle = neighborTable[iNeigh];
		if (iParticle == jParticle) continue;

		xji = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
		if(neighborTablePeriodic[iNeigh]==1){
			xji += parameter.pipeLength;
		}
		else if(neighborTablePeriodic[iNeigh]==2){
			xji -= parameter.pipeLength;
		}
		yji = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);

		distanceIJ_squared = xji*xji + yji*yji;

		if(NumberOfDimensions == 3){
			zji = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
			distanceIJ_squared += zji*zji;
		}

		distanceIJ = sqrt(distanceIJ_squared);

		weight = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForInterface);
		statVolume += weight;

		matrix[XDIM][XDIM] += weight * xji * xji / distanceIJ_squared;
		matrix[XDIM][YDIM] += weight * xji * yji / distanceIJ_squared;
		matrix[YDIM][XDIM] += weight * xji * yji / distanceIJ_squared;
		matrix[YDIM][YDIM] += weight * yji * yji / distanceIJ_squared;
		if (NumberOfDimensions == 3){
			matrix[XDIM][ZDIM] += weight * xji * zji / distanceIJ_squared;
			matrix[ZDIM][XDIM] += weight * xji * zji / distanceIJ_squared;
			matrix[YDIM][ZDIM] += weight * yji * zji / distanceIJ_squared;
			matrix[ZDIM][YDIM] += weight * yji * zji / distanceIJ_squared;
			matrix[ZDIM][ZDIM] += weight * zji * zji / distanceIJ_squared;
		}
	}

	statVolume = NumberOfDimensions / statVolume;
	for(i=0; i < MAX_DIM; i++){
		for(j=0; j < MAX_DIM; j++){
			matrix[i][j] *= statVolume;
		}
	}

	if(NumberOfDimensions == 2) INTERFACE_matrixInverse2x2(matrix,corrMatrixIJ);
	else if(NumberOfDimensions == 3) INTERFACE_matrixInverse3x3(matrix,corrMatrixIJ);

}


void
INTERFACE_matrixInverse2x2b(double inverse[10][10]){

	double determinant;
	double matrix[MAX_DIM][MAX_DIM];
	int i,j;

	for(i=0; i < MAX_DIM; i++){
		for(j=0; j < MAX_DIM; j++){
			matrix[i][j] = inverse[i][j];
		}
	}

	determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

	if (fabs(determinant) < TOO_SMALL_VALUE){
		printf("Warning: Determinant is 0 in MATHE_matrixInverse2x2, setting identity matrix!\n");
		inverse[0][0] = 1;
		inverse[0][1] = 0;
		inverse[1][0] = 0;
		inverse[1][1] = 1;
		return;
	}


	inverse[0][0] = matrix[1][1] / determinant;
	inverse[0][1] = - matrix[0][1] / determinant;
	inverse[1][0] = - matrix[1][0] / determinant;
	inverse[1][1] = matrix[0][0] / determinant;

}

void
INTERFACE_matrixInverse2x2(double matrix[MAX_DIM][MAX_DIM], double inverse[MAX_DIM][MAX_DIM]){

	double determinant;

	determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

	if (fabs(determinant) < TOO_SMALL_VALUE){
		printf("Warning: Determinant is 0 in MATHE_matrixInverse2x2, setting identity matrix!\n");
		inverse[0][0] = 1;
		inverse[0][1] = 0;
		inverse[1][0] = 0;
		inverse[1][1] = 1;
		return;
	}


	inverse[0][0] = matrix[1][1] / determinant;
	inverse[0][1] = - matrix[0][1] / determinant;
	inverse[1][0] = - matrix[1][0] / determinant;
	inverse[1][1] = matrix[0][0] / determinant;

}

void
INTERFACE_matrixInverse3x3(double matrix[MAX_DIM][MAX_DIM], double inverse[MAX_DIM][MAX_DIM]){

	double determinant;
	int i, j;

	determinant = 0;
	for(i = 0; i < 3; i++)
		determinant += (matrix[0][i] * (matrix[1][(i+1)%3] * matrix[2][(i+2)%3] - matrix[1][(i+2)%3] * matrix[2][(i+1)%3]));

	if (fabs(determinant) < TOO_SMALL_VALUE){
		printf("Warning: Determinant is 0 in MATHE_matrixInverse2x2, setting identity matrix!\n");
		for(i = 0; i < 3; i++){
			for(j = 0; j < 3; j++){
				if(i==j) inverse[i][j] = 1;
				else inverse[i][j] = 0;
			}
		}
		return;
	}

	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			inverse[i][j] = ((matrix[(j+1)%3][(i+1)%3] * matrix[(j+2)%3][(i+2)%3]) - (matrix[(j+1)%3][(i+2)%3] * matrix[(j+2)%3][(i+1)%3]))/ determinant;

}



void
INTERFACE_initializeInterfaceFile(void){

	FILE *fp;
	int i,j;

	if (parameter.flagOfInterfaceCalculation == OFF) return;

	if (timer.simulationTime != 0.0)
		if (fp = fopen("Interface.area", "r")){
			fclose(fp);
			return;
		}

	if (parameter.flagOfInterfaceCompare == ON){
		INTERFACE_initializeInterfaceFileCompare();
		return;
	}

	fp = FILE_openFile("Interface.area", "w");

	fprintf(fp, "Time[s]");

	for(i=0;i<parameter.numberOfParticleTypes;i++){

		if (i == GHOST) continue;
		if (i == parameter.dummyWallType) continue;
		if (i == parameter.wallType) continue;

		for (j=0;j<parameter.numberOfParticleTypes;j++){

			if (j == GHOST) continue;
			if (j == parameter.dummyWallType) continue;
			if (j == parameter.wallType) continue;
			if (j == i) continue;

			if (j == parameter.wallType)
				fprintf(fp, "	Interface%dwithWall",i);
			else
				fprintf(fp, "	Interface%dwith%d",i,j);

			if (j == parameter.wallType) continue;
			if (j == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;
		}
		fprintf(fp, "	TotalInterface%d", i);
		fprintf(fp, "	TotalLiquidInterface%d", i);
	}

	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);

}

void
INTERFACE_initializeInterfaceFileCompare(void){

	FILE *fp;

	if (parameter.flagOfInterfaceCalculation == OFF) return;

	fp = FILE_openFile("Interface.area", "w");

	if (parameter.flagOfSingleVortex == ON)
		fprintf(fp, "Time[s]	Conv01Th	Conv10Th	ConvCombTh	Corr01Th	Corr10Th	CorrCombTh	LS01Th	LS10Th	LSCombTh	LS201Th	LS210Th	LS2CombTh");
	else{
		//fprintf(fp, "Time[s]	Uncorr	Corr");
		fprintf(fp, "Time[s]	WallLiquidTh	LiquidWallTh	WallLiquidCombTh");//	LiquidWallComb");
		//fprintf(fp, "Time[s]	Wall0	Wall1	WallAll	Wall0Comb	Wall1Comb	WallAllComb");
	fprintf(fp, "Time[s]	Conv	ConvAbs	Corr	CorrAbs	LSMPSp1	LSMPSp1Abs");//	LSMPSp1	LSMPSp2");
	//fprintf(fp, "Time[s]	Conv	ConvAbs	LSMPSp1	LSMPSp1Abs	LSMPSp2	LSMPSp2Abs");
	//fprintf(fp, "Time[s]	unCor0with1	Cor0with1	Lim0with1	LimCor0with1	unCor1with0	Cor1with0	Lim1with0	LimCor1with0");
	}

	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);

}

void
INTERFACE_writeInterfaceFile(int iType, int jType){

	FILE *fp;
	int   i,j;

	fp = FILE_openFile("Interface.area", "a");

	fprintf(fp, "\n");
	fprintf(fp, "%lf", timer.simulationTime);

	for(i=0;i<parameter.numberOfParticleTypes;i++){

		if (i == GHOST) continue;
		if (i == parameter.dummyWallType) continue;

		for (j=0;j<parameter.numberOfParticleTypes;j++){

			if (j == GHOST) continue;
			if (j == parameter.dummyWallType) continue;
			if (j == i) continue;

			fprintf(fp, "	%lf", parameter.interfaceArea[i][j]);

			if (j == parameter.wallType) continue;
			if (j == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;
		}
		fprintf(fp, "	%lf", parameter.totalInterfaceArea[i]);
		fprintf(fp, "	%lf", parameter.totalLiquidInterfaceArea[i]);
	}

	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);
}


void
INTERFACE_writeInterfaceFileCompare(int iType, int jType){

	FILE *fp;

	fp = FILE_openFile("Interface.area", "a");

	fprintf(fp, "	%lf", parameter.interfaceArea[iType][jType]);

	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);

}
