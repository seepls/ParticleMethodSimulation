/*****************************************************************************

 MPS-SW-MAIN

 Version:          2.0

 Completion date:  February 14, 2008

 Copyright:        Kazuya SHIBATA and Seiichi KOSHIZUKA
 hold the copyright of this code.
 (2013-2015 Added Bubble Model for Stirred Tank Analysis by Shogo KAITO)
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "file.h"
#include "memory.h"
#include "copy.h"
#include "bucket.h"
#include "domain.h"
#include "distance.h"
#include "neigh.h"
#include "density.h"
#include "timer.h"
#include "other.h"
#include "gravity.h"
#include "convection.h"
#include "collision.h"
#include "pressure.h"
#include "gradient.h"
#include "viscosity.h"
#include "sort.h"
#include "object.h"
#include "forcedMotion.h"
#include "init.h"
#include "finalize.h"


#include "bubble.h"
#include "buoyancy.h"
#include "concentration.h"
#include "interface.h"
#include "vortex.h"


int
main( int argumentCount, char **argumentVector ){

    double tic;
	
    tic = TIMER_getWTime();
    INIT_initializeParameters( argumentCount, argumentVector );                                  //initialize parameters
	
    omp_set_num_threads(parameter.threads);

    for( timer.iTimeStep= 0; timer.iTimeStep < timer.finishTimeStep;  timer.iTimeStep++ ){


        TIMER_setDtAutomatically();                                                              //set Dt
		//timer.dt = timer.dt_initial; // Use only for easy simulations
		TIMER_putTimeForwardByDt();                                                              //put timer forward by Dt

        TIMER_displayStateOfTimeStep_atAppropriateTime();                                        //display state of Timestep

        GRAVITY_calculateGravity();                                                              //calculate Gravity

        VISCOSITY_calculateViscosity();                                                          //calculate Viscosity

        VORTEX_setVelocityOfVortex();

		CONVECTION_moveParticles( particle.position, particle.velocity );                        //move Particles by convection
		
		VORTEX_redistributeParticles();
		
        NEIGH_setNeighborTable(particle.position);                                               //set Neighbor Table

		if(parameter.flagOfForcedMotionOfRigidBody == ON ){                                      //main Function of Forced Motion
            FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( &forcedMotionOfRigidBody );
        }

        COLLISION_calculateCollisionBetweenParticles();                                          //calculate Collision between Particles

        OTHER_checkThatParticlesAreNotAllGhost();                                                //check that Particles are not all ghost

        PRESSURE_calculatePressure();                                                            //calculate Pressure

        GRADIENT_correctParticleVelocityAndPositionUsingPressureGradient();                      //correct velocity and position

        if(parameter.flagOfForcedMotionOfRigidBody == ON ){                                      //main Function of Forced Motion
            FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( &forcedMotionOfRigidBody );
        }

        COPY_updateParticleProperty();                                                           //update Particle Properties

        if(parameter.flagOfBubbleCalculation == ON){                                             //calculate Bubbles
            Bubble_calculateBubble();
        }

		if (parameter.flagOfConcentrationCalculation == ON){                                     //calculate Concentration
			Concentration_calculateConcentration();
		}

        COLLISION_calculateCollisionBetweenParticles();                                          //calculate Collision between Particles

        NEIGH_setNeighborTable(particle.position);                                               //set Neighbor Table
		
        INTERFACE_CalculateInterfaceArea();

        FILE_writeCalculationResultInFile();                                                     //write Result in File

        if( YES == TIMER_checkWhetherItIsTimeToFinishProgram() ) break;                          //check finishing time

    }

    printf("Time for Execution: %f s\n",TIMER_getWTime() - tic);

    FINALIZE_finalizeProgram();                                                                  //finalize program

    return(EXIT_SUCCESS);

}



