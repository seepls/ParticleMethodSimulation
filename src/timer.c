#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "maxmin.h"
#include "other.h"
#include "timer.h"


void
TIMER_initializeTimeFunctions( void ){


  if( timer.typeOfOutputInterval == SIMULATION_TIME ){

	if( timer.intervalTimeOfWritingProfFile_simulationTime > 0.0){
	  timer.initialSequantialNumberOfProfFile = (int)(timer.simulationTime / timer.intervalTimeOfWritingProfFile_simulationTime);

	  fprintf(FpForLog,"timer.initialSequantialNumberOfProfFile= %d\n",timer.initialSequantialNumberOfProfFile);

	}else{
	  OTHER_endProgram("timer.intervalTimeOfWritingProfFile_simulationTime is wrong.");
	}

  }else{

	timer.initialSequantialNumberOfProfFile = 0;

  }

  timer.iProfFile = timer.initialSequantialNumberOfProfFile;

  timer.iVtkFile = timer.initialSequantialNumberOfProfFile;


  timer.iPressureFile = 0;

  timer.maxDt = timer.maxDt_ratio * timer.dt_initial;

  timer.minDt = timer.minDt_ratio * timer.dt_initial;

  timer.timeToWriteProfFile_simulationTime = timer.simulationTime + timer.intervalTimeOfWritingProfFile_simulationTime;

  timer.timeToWriteProfFile_timeStep       = 0;

  /*
  timer.timeToWriteProfFile_timeStep       = 0 + timer.intervalTimeOfWritingProfFile_timeStep;
  */


  timer.timeToDisplayStateOfTimeStep       = 0.0;



}





void
TIMER_setDtAutomatically( void ){

  double dt_previous;
  double upperLimitOfDt;
  double smallerDt;

  double dt_byCourantNumber;
  double dt_byDiffusionNumber;


  if(timer.iTimeStep == 0){
    timer.dt = timer.dt_initial;
  }


  dt_previous        = timer.dt;
  dt_byCourantNumber = TIMER_calculateDtWithCourantNumberInMind();


  if( parameter.flagOfViscosityCalculation == ON && parameter.flagOfHighViscosityCalculation != ON){
    dt_byDiffusionNumber = TIMER_calculateDtWithDiffusionNumberInMind();
  }else{
    dt_byDiffusionNumber = dt_byCourantNumber;
  }


  upperLimitOfDt = dt_previous * timer.upperLimitOfChangeRateOfDt;

  smallerDt      = MAXMIN_returnSmallerValue( dt_byCourantNumber, dt_byDiffusionNumber );
  smallerDt      = MAXMIN_returnSmallerValue( smallerDt,          upperLimitOfDt       );
  smallerDt      = MAXMIN_returnSmallerValue( smallerDt,          timer.maxDt          );

  if( timer.iTimeStep == 0){
    smallerDt = MAXMIN_returnSmallerValue(smallerDt, timer.dt_initial);
  }

  timer.dt = smallerDt;

  TIMER_checkThatDtIsNotTooSmall();

}



double
TIMER_calculateDtWithCourantNumberInMind( void ){

  double dt;

  MAXMIN_calculateMaxAbsoluteValueOfVelocity();

  if( particle.maxAbsoluteValueOfVelocity > 0.0){

    dt = (parameter.courantNumber * particle.averageDistance)/particle.maxAbsoluteValueOfVelocity;

  }else{

    dt = timer.dt_initial;

  }

  return dt;

}




double
TIMER_calculateDtWithDiffusionNumberInMind( void ){

  double dt;

  dt = parameter.diffusionNumber * (particle.averageDistance * particle.averageDistance) / (2.0 * physicalProperty.kinematicViscosity);

  return dt;

}


void
TIMER_checkThatDtIsNotTooSmall( void ){

  char errorMessage[256];

  if( timer.dt < timer.minDt ) {
    sprintf(errorMessage,"[in TIMER_setDtAutomatically]  dt is too small.\n dt = %lf  maxAbsoluteValueOfVelocity = %lf",timer.dt, particle.maxAbsoluteValueOfVelocity);

    OTHER_endProgram(errorMessage);
  }

}



void
TIMER_putTimeForwardByDt( void ){

  timer.simulationTime += timer.dt;

}



int
TIMER_checkWhetherItIsTimeToWriteProfFile( void ){

  int answer = NO;

  if( timer.typeOfOutputInterval == SIMULATION_TIME ){

    if( timer.simulationTime >= timer.timeToWriteProfFile_simulationTime ){
      timer.timeToWriteProfFile_simulationTime += timer.intervalTimeOfWritingProfFile_simulationTime;
      answer = YES;
    }else{
      answer = NO;
    }

  }else if(timer.typeOfOutputInterval == TIME_STEP ){

    if( timer.iTimeStep >= timer.timeToWriteProfFile_timeStep ){
      timer.timeToWriteProfFile_timeStep += timer.intervalTimeOfWritingProfFile_timeStep;
      answer = YES;
    }else{
      answer = NO;
    }

  }else{

    answer = NO;

  }

  return(answer);

}




int
TIMER_checkWhetherItIsTimeToFinishProgram( void ){

  if( timer.simulationTime > timer.finishTime){
    return(YES);
  }

  return(NO);

}



void
TIMER_displayStateOfTimeStep_atAppropriateTime( void ){

  timer.flagOfDisplayingStateOfTimeStep = OFF;

  if( YES == TIMER_checkWhetherItIsTimeToDisplayStateOfTimeStep() ){
    if (parameter.flagOfHideOutput == OFF){
        TIMER_displayStateOfTimeStep( stderr );

        timer.flagOfDisplayingStateOfTimeStep = ON;
    }
  }

  TIMER_displayStateOfTimeStep( FpForLog );


}



int
TIMER_checkWhetherItIsTimeToDisplayStateOfTimeStep( void ){

  int answer = NO;


  if( (timer.calculationTime >= timer.timeToDisplayStateOfTimeStep)|| (timer.iTimeStep <= 100) ){

	timer.timeToDisplayStateOfTimeStep = timer.calculationTime + timer.intervalTimeOfDisplayingStateOfTimeStep;
	answer = YES;

  }else{
	answer = NO;
  }

  return(answer);

}


void
TIMER_displayStateOfTimeStep( FILE *fp ){

  TIMER_calculateCalculationTime();
  TIMER_displaySimulationTime( fp );
  TIMER_displayCalculationTime( fp );

}




void
TIMER_displaySimulationTime( FILE *fp ){

  fprintf(fp,"\n");
  fprintf(fp,"--------------------------------------------------------------------------------\n");
  fprintf(fp,">> iTimeStep= %d,  numberOfParticles= %d,   dt= %e[sec],  simulationTime= %lf[sec]\n"
	  ,timer.iTimeStep, particle.totalNumber, timer.dt, timer.simulationTime);

}



void
TIMER_displayCalculationTime( FILE *fp ){

  //if( parameter.flagOfHideOutput == ON) return;

  fprintf(fp,"\n");

  fprintf(fp,"calculationTime ---- total:");

  if( timer.hours > 0 ){

    fprintf(fp,"   %d[hour]   %d[minute]  %.3lf[second]\n",timer.hours, timer.minutes, timer.seconds );

  }else if( timer.minutes > 0 ){

    fprintf(fp,"   %d[minute]  %.3lf[second]\n",timer.minutes, timer.seconds );

  }else{

    fprintf(fp," %10lf [ second      ]\n", timer.seconds );
  }


  if( timer.iTimeStep >= 1){
    fprintf(fp,"                ---- 1step: %10lf [ second/step ]\n",timer.calculationTimePerTimeStep );
  }

  TIMER_displayRealTime(fp);


}

void
TIMER_displayRealTime( FILE *fp ){
	if(parameter.threads == 1) return;

  fprintf(fp,"\n");

  fprintf(fp,"realTime ---- total:");

  if( timer.realHours > 0 ){

    fprintf(fp,"   %d[hour]   %d[minute]  %.3lf[second]\n",timer.realHours, timer.realMinutes, timer.realSeconds );

  }else if( timer.realMinutes > 0 ){

    fprintf(fp,"   %d[minute]  %.3lf[second]\n",timer.realMinutes, timer.realSeconds );

  }else{

    fprintf(fp," %10lf [ second      ]\n", timer.realSeconds );
  }


  if( timer.iTimeStep >= 1){
    fprintf(fp,"         ---- 1step: %10lf [ second/step ]\n",timer.realTimePerTimeStep );
  }


}



void
TIMER_calculateCalculationTime( void ){

  timer.clockAtCurrentTime = clock();
  timer.realTimeCurrent = TIMER_getWTime();

  timer.calculationTime = (double)(timer.clockAtCurrentTime - timer.clockAtStartingTime)/(double)(CLOCKS_PER_SEC);
  timer.realTime = timer.realTimeCurrent - timer.realTimeStart;


  if( timer.iTimeStep >= 1){
	timer.calculationTimePerTimeStep = (double)(timer.clockAtCurrentTime - timer.clockAtBeginningOfMainLoop)/(double)(CLOCKS_PER_SEC);
	timer.calculationTimePerTimeStep /= (double)( timer.iTimeStep);
    timer.realTimePerTimeStep = timer.realTimeCurrent - timer.realTimeLoop;
    timer.realTimePerTimeStep /= (double)(timer.iTimeStep);
  }else{
	timer.calculationTimePerTimeStep = (double)(timer.clockAtCurrentTime - timer.clockAtBeginningOfMainLoop)/(double)(CLOCKS_PER_SEC);
	timer.realTimePerTimeStep = timer.realTimeCurrent - timer.realTimeLoop;
  }


  timer.hours   = (int)(timer.calculationTime/3600.0);
  timer.minutes = (int)( (timer.calculationTime - 3600.0 * timer.hours)/60.0 );
  timer.seconds =  fmod(timer.calculationTime,60.0);

  timer.realHours = (int)(timer.realTime/3600.0);
  timer.realMinutes = (int)((timer.realTime - 3600 * timer.realHours)/60.0);
  timer.realSeconds = fmod(timer.realTime,60.0);

}


double
TIMER_getWTime(void){
	//struct timeval t;
#if defined(_OPENMP)
	return omp_get_wtime();
#else
	//struct timeval t;

	//gettimeofday(&t, NULL);
	//return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
#endif
}

