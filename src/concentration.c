#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "weight.h"
#include "extern.h"
#include "concentration.h"
#include "buoyancy.h"
#include "distance.h"
#include "file.h"
#include "gravity.h"
#include "neigh.h"

void
Concentration_calculateConcentration(void){

	if (parameter.flagOfConcentrationCalculation == OFF) return;
	              
	particle.concentration_previous = particle.concentration;

	Concentration_calculateConcentrationOfLiquidParticles();

	Concentration_checkSumOfConcentrations();

	Concentration_calculateMixingIndex();

	Concentration_writeMixingIndexFile();

}


void
Concentration_calculateConcentrationOfLiquidParticles(void){


	
	int iParticle, jParticle;
	int iDim;
	int iNeigh;

	double sigma[MAX_DIM];
	double distanceIJ;
	double weightIJ;

	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	double DiffusionCoefficient;
	double radiusOfLaplacianForDiffusion;
	double lambdaTimesNZero;
	double dt;
	
	if (parameter.flagOfConcentrationCalculation == OFF) return;

	DiffusionCoefficient = parameter.DiffusionCoefficient;
	radiusOfLaplacianForDiffusion = parameter.radiusOfLaplacianForViscosity;
	lambdaTimesNZero = parameter.lambdaOfLaplacianForDiffusion * parameter.nZeroOfLaplacianForDiffusion;
	dt = timer.dt;

	NEIGH_selectNeighborTable(&numberOfNeighborParticles
		, &neighborTable
		, parameter.radiusOfLaplacianForDiffusion_ratio
		,&neighborTablePeriodic
		);


	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){

		if (particle.type[iParticle] == GHOST) continue;
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
		if (particle.type[iParticle] == parameter.wallType) continue;
		if (particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;


		for (iDim = 0; iDim < NumberOfDimensions; iDim++){
			sigma[iDim] = 0.0;
		}

		for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			if (particle.type[jParticle] == GHOST)continue;
			if (particle.type[jParticle] == parameter.dummyWallType) continue;
			if (particle.type[jParticle] == parameter.wallType) continue;
			if (particle.type[jParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;

			//distanceIJ = DISTANCE_calculateDistanceBetweenParticles(particle.position, iParticle, jParticle);
			distanceIJ = DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position, iParticle, jParticle,neighborTablePeriodic[iParticle][iNeigh]);
			weightIJ = WEIGHT_calculateWeightFunction(distanceIJ, radiusOfLaplacianForDiffusion);

			for (iDim = 0; iDim < NumberOfDimensions; iDim++){
				sigma[iDim] += (particle.concentration_previous[iDim][jParticle] - particle.concentration_previous[iDim][iParticle]) *  weightIJ;
			}

		}
		for (iDim = 0; iDim < NumberOfDimensions; iDim++){

			sigma[iDim] *= (2.0 * NumberOfDimensions) / lambdaTimesNZero;
			particle.concentration[iDim][iParticle] += DiffusionCoefficient * sigma[iDim] * dt;
		}
	}
}

void
Concentration_checkSumOfConcentrations(void){

	int iParticle;
	int iType;
	double sum;
	double difference;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){

		if (particle.type[iParticle] == GHOST) continue;
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
		if (particle.type[iParticle] == parameter.wallType) continue;
		if (particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;

		sum = particle.concentration[0][iParticle] + particle.concentration[1][iParticle] + particle.concentration[2][iParticle];

		if (sum == 1.0) continue;

		difference = 1.0 - sum;

		for (iType = 0; iType < 3; iType++){

		particle.concentration[iType][iParticle] += difference * (particle.concentration[iType][iParticle] / sum);

		}
	}
}

void
Concentration_averageConcentration(void){

	int iType;

	for (iType = 0; iType < 3; iType++){
		parameter.averageConcentration[iType] = (double)parameter.numberOfLiquidParticles[iType] / (double)parameter.numberOfFluidParticles;
	}

	for (iType = 0; iType < 3; iType++){
		parameter.maxConcentrationDifference[iType] = parameter.numberOfLiquidParticles[iType] * (1 - parameter.averageConcentration[iType])
			+ (parameter.numberOfFluidParticles - parameter.numberOfLiquidParticles[iType]) * parameter.averageConcentration[iType];
		parameter.maxConcentrationDifference[iType] = parameter.maxConcentrationDifference[iType] / parameter.numberOfFluidParticles;
	}
}

void
Concentration_calculateMixingIndex(void){

	int iParticle;
	int iType;
	double sum[3] = { 0, 0, 0 };

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){
		if (particle.type[iParticle] == GHOST) continue;
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
		if (particle.type[iParticle] == parameter.wallType) continue;
		if (particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;

		for (iType = 0; iType < 3; iType++){
			sum[iType] += fabs(parameter.averageConcentration[iType] - particle.concentration[iType][iParticle]);
		}
	}
	for (iType = 0; iType < 3; iType++){
		if (parameter.maxConcentrationDifference[iType] != 0) parameter.mixingIndex[iType] = 1 - sum[iType] / (parameter.numberOfFluidParticles * parameter.maxConcentrationDifference[iType]);
		else parameter.mixingIndex[iType] = 0;
	}
}


void
Concentration_writeMixingIndexFile(void){

	FILE *fp;
	int   iType;

	if (parameter.flagOfConcentrationCalculation == OFF) return;

	fp = FILE_openFile(parameter.nameOfMixingIndexFile, "a");

	fprintf(fp, "\n");
	fprintf(fp, "%lf", timer.simulationTime);
	for (iType = 0; iType < 3; iType++){
		fprintf(fp, "	%lf", parameter.mixingIndex[iType]);
	}



	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);

}

void
Concentration_initializeMixingIndexFile(void){

	FILE *fp;
	int   iType;

	if (parameter.flagOfConcentrationCalculation == OFF) return;

	fp = FILE_openFile(parameter.nameOfMixingIndexFile, "w");

	fprintf(fp, "Time	Liquid1	Liquid2	Liquid3\n");

	fprintf(fp, "%lf", timer.simulationTime);
	for (iType = 0; iType < 3; iType++){
		fprintf(fp, "	%lf", parameter.mixingIndex[iType]);
	}


	fflush(fp);
	FILE_closeFile(fp, parameter.nameOfMixingIndexFile);

}
