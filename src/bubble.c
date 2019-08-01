#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "weight.h"
#include "extern.h"
#include "bubble.h"
#include "buoyancy.h"
#include "distance.h"
#include "file.h"
#include "gravity.h"
#include "neigh.h"

void
Bubble_calculateBubble(void){

    Bubble_calculateBubbleRising();

    Bubble_surfaceJudgement();

    Buoyancy_calculateBuyancy();

}

double
Bubble_calculateCosine(double distanceIJ,int iParticle,int jParticle){

    double cosine;

    cosine=(particle.position[YDIM][jParticle]-particle.position[YDIM][iParticle])/fabs(distanceIJ);

    return(cosine);
}

double
Bubble_setInfluenceRadius(void){

    double influenceRadius;

    if(NumberOfDimensions == 2){

        influenceRadius=4.1*particle.averageDistance;

    }else{

        influenceRadius=3.1*particle.averageDistance;

    }

    return(influenceRadius);
}

double
Bubble_calculateVolume(double diameter){

    double volume;

    volume=(4.0/3.0)*M_PI*pow((diameter/2.0),3.0);

    return(volume);
}

void
Bubble_setBetaZero(void){

    int iParticle;
    double influenceRadius;
    double distanceIJ;
    double cosine;
    double Beta;

    Beta=0.0;

    influenceRadius=Bubble_setInfluenceRadius();

    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){

        if(particle.type[iParticle]==parameter.wallType || particle.type[iParticle]== parameter.dummyWallType)continue;

        distanceIJ=DISTANCE_calculateDistanceBetweenParticles(particle.position,iParticle,parameter.numperOfParticleForCalculatingBeta);

        if(distanceIJ>influenceRadius||iParticle==parameter.numperOfParticleForCalculatingBeta)continue;

        cosine=(particle.position[YDIM][iParticle]-particle.position[YDIM][parameter.numperOfParticleForCalculatingBeta])/distanceIJ;

        if(cosine<=0.0)continue;

        Beta+=(particle.position[YDIM][iParticle]-particle.position[YDIM][parameter.numperOfParticleForCalculatingBeta])*cosine*WEIGHT_calculateWeightFunction(distanceIJ,influenceRadius);

    }

    parameter.betaZeroOfParticles=Beta;
}

double
Bubble_calculateDiameter(int iParticle){

    double volume;
    double diameter;
    double voidrate;


    if(particle.numberOfBubbles[iParticle]>0.0){

        volume=Bubble_calculateVolume(particle.averageDistance);
        voidrate=(particle.moleOfBubbles[iParticle]*physicalProperty.gasConstant*physicalProperty.temperature)/((particle.pressure[iParticle]+physicalProperty.headPressure)*volume);
        diameter=particle.averageDistance*pow((voidrate/particle.numberOfBubbles[iParticle]),1.0/3.0);

    }else{

        diameter=0.0;

    }

    return(diameter);
}

double
Bubble_calculateRisingVelocity(double diameter){

    double risingVelocity;
    double viscosityCoefficientofFluid;
    double Re;

    viscosityCoefficientofFluid=physicalProperty.massDensity[0]*physicalProperty.kinematicViscosity;

    risingVelocity=(pow(diameter,2)*(physicalProperty.massDensityOfBubble-physicalProperty.massDensity[0]))*physicalProperty.gravity[YDIM]/(18.0*viscosityCoefficientofFluid);

    Re=(risingVelocity*diameter)/physicalProperty.kinematicViscosity;

    if(2<=Re){

        risingVelocity=pow(((4.0/225.0)*pow(((physicalProperty.massDensityOfBubble-physicalProperty.massDensity[0])*physicalProperty.gravity[YDIM]),2.0)/(physicalProperty.massDensity[0]*viscosityCoefficientofFluid)),1.0/3.0)*diameter;
        Re=(risingVelocity*diameter)/physicalProperty.kinematicViscosity;

        if(500<=Re){

            risingVelocity=pow((4.0/(3*0.44))*(physicalProperty.massDensityOfBubble-physicalProperty.massDensity[0])*physicalProperty.gravity[YDIM]*diameter/physicalProperty.massDensity[0],1.0/2.0);
            Re=(risingVelocity*diameter)/physicalProperty.kinematicViscosity;
        }
    }

    return(risingVelocity);
}


void
Bubble_calculateBubbleRising(void){

    int iParticle,jParticle;
    double risingVelocity;
    double distanceIJ;
    double cosine;
    double influenceRadius;
    double beta;
    double aIJ;
    double bubbleDiameter;

    int    iNeigh;
    int    *numberOfNeighborParticles;
    int    **neighborTable;
	int    **neighborTablePeriodic;

    NEIGH_selectNeighborTable(  &numberOfNeighborParticles
                              ,&neighborTable
                              ,parameter.radiusOfLaplacianForPressure_ratio
							  ,&neighborTablePeriodic
                              );

    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){

        if(particle.type[iParticle]==parameter.wallType || particle.type[iParticle]== parameter.dummyWallType ||particle.type[iParticle]== parameter.typeNumberOfRigidParticle_forForcedMotion)continue;


        bubbleDiameter=Bubble_calculateDiameter(iParticle);
        risingVelocity=Bubble_calculateRisingVelocity(bubbleDiameter);

        if(risingVelocity<=0.0)continue;

        beta=parameter.betaZeroOfParticles/risingVelocity;
        influenceRadius=Bubble_setInfluenceRadius();


        for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){

            jParticle = neighborTable[iParticle][iNeigh];

            if(particle.type[jParticle]==parameter.wallType || particle.type[jParticle]== parameter.dummyWallType ||particle.type[jParticle]== parameter.typeNumberOfRigidParticle_forForcedMotion)continue;

            //distanceIJ=DISTANCE_calculateDistanceBetweenParticles(particle.position,iParticle,jParticle);
			distanceIJ=DISTANCE_calculateDistanceBetweenParticlesPeriodic(particle.position,iParticle,jParticle,neighborTablePeriodic[iParticle][iNeigh]);
            cosine=Bubble_calculateCosine(distanceIJ,iParticle,jParticle);

            if(cosine<=0.0||jParticle==iParticle)continue;

            aIJ=(timer.dt/beta)*WEIGHT_calculateWeightFunction(distanceIJ,influenceRadius)*cosine;
            particle.numberOfBubbles[jParticle]+= aIJ*particle.numberOfBubbles_previous[iParticle];
            particle.numberOfBubbles[iParticle]-= aIJ*particle.numberOfBubbles_previous[iParticle];
            particle.moleOfBubbles[jParticle]+= aIJ*particle.moleOfBubbles_previous[iParticle];
            particle.moleOfBubbles[iParticle]-= aIJ*particle.moleOfBubbles_previous[iParticle];

        }
    }
}

void
Bubble_surfaceJudgement(void){

    int iParticle;
    double timeConstant;
    double risingVelocity;
    double diameter;


    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){

        if(particle.type[iParticle]==parameter.wallType || particle.type[iParticle]== parameter.dummyWallType)continue;

        if(particle.particleNumberDensity[iParticle]< 0.97*parameter.nZeroOfParticleNumberDensity){

            diameter=Bubble_calculateDiameter(iParticle);
            risingVelocity=Bubble_calculateRisingVelocity(diameter);
            timeConstant=(particle.averageDistance)/risingVelocity;

            particle.numberOfBubbles[iParticle]-=(timer.dt/timeConstant)*particle.numberOfBubbles[iParticle];
            particle.moleOfBubbles[iParticle]-=(timer.dt/timeConstant)*particle.moleOfBubbles[iParticle];

            if(particle.numberOfBubbles[iParticle]<=0.0 || particle.moleOfBubbles[iParticle]<=0.0){

               particle.numberOfBubbles[iParticle]=0.0;
               particle.moleOfBubbles[iParticle]=0.0;

            }
        }
    }
}


