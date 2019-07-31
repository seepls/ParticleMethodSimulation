#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "struct.h"
#include "copy.h"
#include "neigh.h"
#include "distance.h"
#include "extern.h"



void
VORTEX_setVelocityOfVortex(void){

	int iParticle;

	if (parameter.flagOfSingleVortex == OFF) return;

#pragma omp parallel for
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){
		particle.velocity[XDIM][iParticle]=sin(M_PI*particle.position[0][iParticle])*cos(M_PI*particle.position[1][iParticle]);
		particle.velocity[YDIM][iParticle]=-cos(M_PI*particle.position[0][iParticle])*sin(M_PI*particle.position[1][iParticle]);
	}

	/*if( timer.simulationTime > timer.finishTime/2){
		for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){
			particle.velocity[0][iParticle] *=-1;
			particle.velocity[1][iParticle] *=-1;
		}
	}*/
}

void
VORTEX_redistributeParticles(void){

	int iParticle, jParticle;
	int iNeigh;

	double xij;
	double yij;
	//double nxij;
	//double nyij;
	double sumx;
	double sumy;

	double distanceIJ;
	double distanceIJ_squared;

	int    *numberOfNeighborParticles;
	int    **neighborTable;
	int    **neighborTablePeriodic;

	double drix;
	double driy;
	double Ri_squared;

	double mirrorMin[MAX_DIM];
	double mirrorMax[MAX_DIM];

	if (parameter.flagOfSingleVortex == OFF) return;

	mirrorMin[XDIM] = 0;// + particle.averageDistance/2;// + parameter.radiusOfGradient;
	mirrorMax[XDIM] = 1;// - particle.averageDistance/2;// - parameter.radiusOfGradient;
	mirrorMin[YDIM] = 0;// + particle.averageDistance/2;// + parameter.radiusOfGradient;
	mirrorMax[YDIM] = 1;// - particle.averageDistance/2;// - parameter.radiusOfGradient;

	COPY_copy2dimDoubleArray( 3, particle.totalNumber, particle.position_previous,      particle.position);

	NEIGH_selectNeighborTable(   &numberOfNeighborParticles
			       ,&neighborTable
			       ,parameter.maxRadius
				   ,&neighborTablePeriodic
			       );
#pragma omp parallel for private(jParticle, iNeigh, sumx, sumy, Ri_squared, drix, driy, xij, yij, distanceIJ, distanceIJ_squared)
	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
		if (particle.type[iParticle] == GHOST) continue;
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
		if (particle.type[iParticle] == parameter.wallType) continue;
		if (particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion && parameter.flagOfForcedMotionOfRigidBody == ON) continue;

		sumx = 0;
		sumy = 0;
		Ri_squared = 0;

		for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			xij = (particle.position_previous[XDIM][jParticle] - particle.position_previous[XDIM][iParticle]);
			yij = (particle.position_previous[YDIM][jParticle] - particle.position_previous[YDIM][iParticle]);

			distanceIJ_squared = xij*xij + yij*yij;
			distanceIJ = sqrt(distanceIJ_squared);

			sumx += (xij/distanceIJ) / distanceIJ_squared;
			sumy += (yij/distanceIJ) / distanceIJ_squared;
			Ri_squared += distanceIJ;

		}

		if ((particle.position_previous[XDIM][iParticle] <= mirrorMin[XDIM] + parameter.maxRadius)||(particle.position_previous[XDIM][iParticle] >= mirrorMax[XDIM] - parameter.maxRadius)){
			for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle] + 1; iNeigh++){
				if (iNeigh == numberOfNeighborParticles[iParticle])
					jParticle = iParticle;
				else
					jParticle = neighborTable[iParticle][iNeigh];

				if (particle.position_previous[XDIM][iParticle] <= mirrorMin[XDIM] + parameter.maxRadius)
					xij = (2 * mirrorMin[XDIM] - particle.position_previous[XDIM][jParticle] - particle.position_previous[XDIM][iParticle]);
				else if (particle.position_previous[XDIM][iParticle] >= mirrorMax[XDIM] - parameter.maxRadius)
					xij = (2 * mirrorMax[XDIM] - particle.position_previous[XDIM][jParticle] - particle.position_previous[XDIM][iParticle]);
				else continue;
				yij = (particle.position_previous[YDIM][jParticle] - particle.position_previous[YDIM][iParticle]);

				distanceIJ_squared = xij*xij + yij*yij;
				distanceIJ = sqrt(distanceIJ_squared);

				if (distanceIJ >= parameter.maxRadius) continue;
				sumx += (xij/distanceIJ) / distanceIJ_squared;
				sumy += (yij/distanceIJ) / distanceIJ_squared;
				Ri_squared += distanceIJ;
			}
		}
		if ((particle.position_previous[YDIM][iParticle] <= mirrorMin[YDIM] + parameter.maxRadius)||(particle.position_previous[YDIM][iParticle] >= mirrorMax[YDIM] - parameter.maxRadius)){
			for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle] + 1; iNeigh++){
				if (iNeigh == numberOfNeighborParticles[iParticle])
					jParticle = iParticle;
				else
					jParticle = neighborTable[iParticle][iNeigh];

				xij = (particle.position_previous[XDIM][jParticle] - particle.position_previous[XDIM][iParticle]);
				if (particle.position_previous[YDIM][iParticle] <= mirrorMin[YDIM] + parameter.maxRadius)
					yij = (2 * mirrorMin[YDIM] - particle.position_previous[YDIM][jParticle] - particle.position_previous[YDIM][iParticle]);
				else if (particle.position_previous[YDIM][iParticle] >= mirrorMax[YDIM] - parameter.maxRadius)
					yij = (2 * mirrorMax[YDIM] - particle.position_previous[YDIM][jParticle] - particle.position_previous[YDIM][iParticle]);
				else continue;

				distanceIJ_squared = xij*xij + yij*yij;
				distanceIJ = sqrt(distanceIJ_squared);

				if (distanceIJ >= parameter.maxRadius) continue;
				sumx += (xij/distanceIJ) / distanceIJ_squared;
				sumy += (yij/distanceIJ) / distanceIJ_squared;
				Ri_squared += distanceIJ;
			}
		}
		if (((particle.position_previous[XDIM][iParticle] <= mirrorMin[XDIM] + parameter.maxRadius)||(particle.position_previous[XDIM][iParticle] >= mirrorMax[XDIM] - parameter.maxRadius))&&
			((particle.position_previous[YDIM][iParticle] <= mirrorMin[YDIM] + parameter.maxRadius)||(particle.position_previous[YDIM][iParticle] >= mirrorMax[YDIM] - parameter.maxRadius))){
			for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle] + 1; iNeigh++){
				if (iNeigh == numberOfNeighborParticles[iParticle])
					jParticle = iParticle;
				else
					jParticle = neighborTable[iParticle][iNeigh];

				if (particle.position_previous[XDIM][iParticle] <= mirrorMin[XDIM] + parameter.maxRadius)
					xij = (2 * mirrorMin[XDIM] - particle.position_previous[XDIM][jParticle] - particle.position_previous[XDIM][iParticle]);
				else if (particle.position_previous[XDIM][iParticle] >= mirrorMax[XDIM] - parameter.maxRadius)
					xij = (2 * mirrorMax[XDIM] - particle.position_previous[XDIM][jParticle] - particle.position_previous[XDIM][iParticle]);
				else continue;
				if (particle.position_previous[YDIM][iParticle] <= mirrorMin[YDIM] + parameter.maxRadius)
					yij = (2 * mirrorMin[YDIM] - particle.position_previous[YDIM][jParticle] - particle.position_previous[YDIM][iParticle]);
				else if (particle.position_previous[YDIM][iParticle] >= mirrorMax[YDIM] - parameter.maxRadius)
					yij = (2 * mirrorMax[YDIM] - particle.position_previous[YDIM][jParticle] - particle.position_previous[YDIM][iParticle]);
				else continue;

				distanceIJ_squared = xij*xij + yij*yij;
				distanceIJ = sqrt(distanceIJ_squared);

				if (distanceIJ >= parameter.maxRadius) continue;
				sumx += (xij/distanceIJ) / distanceIJ_squared;
				sumy += (yij/distanceIJ) / distanceIJ_squared;
				Ri_squared += distanceIJ;
			}
		}

		Ri_squared = Ri_squared * Ri_squared;

		drix = 0 - sumx * Ri_squared * timer.dt * parameter.alphaZero;
		driy = 0 - sumy * Ri_squared * timer.dt * parameter.alphaZero;

		particle.position[XDIM][iParticle] += drix;
		particle.position[YDIM][iParticle] += driy;

	}

}
