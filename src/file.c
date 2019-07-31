#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "memory.h"
#include "timer.h"
#include "other.h"
#include "file.h"
#include "bubble.h"
#include "mathe.h"


FILE*
FILE_openFile(char *filename, char *mode){

  FILE *fp;
  char errorMessage[256];

  if ((fp=fopen(filename,mode))==NULL) {
    sprintf(errorMessage,"File-open of \"%s\" failed\n\n",filename);
    OTHER_endProgram(errorMessage);
  }

  return(fp);

}



void
FILE_closeFile( FILE *fp, char *fileName){

  char errorMessage[256];

  if ( fp== NULL) {

    sprintf(errorMessage,"The file,\"%s\", is not able to be closed because the file pointer is Null.\n\n",fileName);
    OTHER_endProgram(errorMessage);

  }else{

    fclose( fp );

  }

}





void
FILE_readInputFiles( int argumentCount,char **argumentVector ){


  FILE_setNameOfInputFile( argumentCount, argumentVector );

  FILE_setNameOfOutputFile( argumentCount, argumentVector );

  FILE_openLogFile();

  FILE_displayArgumentOfCommandLine( argumentCount, argumentVector );

  FILE_readDataFile();

  FILE_setNameOfInputFile( argumentCount, argumentVector );

  FILE_displayReadDataFile();

  FILE_readGridFile();

  FILE_readDesignationFileForWritingPressureFile();

  /*
  FILE_countNumberOfParticlesEachType();

  FILE_displayNumberOfParticles();
  */


}




void
FILE_openLogFile( void ){

  FpForLog = FILE_openFile( parameter.nameOfLogFile, "w");

}




void
FILE_setNameOfInputFile( int argumentCount,char **argumentVector ){

  int i;

  FILE_setDefaultNameOfInputFile();

  for(i = 1; i < argumentCount; i++ ) {
        if (strcmp( argumentVector[i], "-inData" ) == 0 ) {
            strcpy(parameter.nameOfDataFile,argumentVector[i+1]);
            i++;
        }
        if (strcmp(argumentVector[i], "-inGrid") == 0){
            strcpy(parameter.nameOfGridFile,argumentVector[i+1]);
            i++;
        }
        if (strcmp(argumentVector[i], "-inBubble") == 0){
            strcpy(parameter.nameOfBubbleInputFile,argumentVector[i+1]);
            i++;
        }
        if (strcmp(argumentVector[i], "-inConc") == 0){
            strcpy(parameter.nameOfConcentrationInputFile,argumentVector[i+1]);
            i++;
        }
    }
/*
  if(argumentCount <= 1) return;

  strcpy(parameter.nameOfDataFile, argumentVector[1]);

  if(argumentCount <= 2) return;

  strcpy(parameter.nameOfGridFile, argumentVector[2]);


  if(argumentCount <= 6) return;

  strcpy(parameter.nameOfBubbleInputFile, argumentVector[6]);

  if (argumentCount <= 9) return;

  strcpy(parameter.nameOfConcentrationInputFile, argumentVector[9]);
  */

}



void
FILE_setNameOfOutputFile( int argumentCount,char **argumentVector ){

  int i;

  FILE_setDefaultNameOfOutputFile();

  for(i = 1; i < argumentCount; i++ ) {
        if (strcmp( argumentVector[i], "-outProf" ) == 0 ) {
            if( parameter.flagOfDivisionOfProfFile == ON){
                strcpy(parameter.nameOfDividedProfFile, argumentVector[i+1]);
            }else{
                strcpy(parameter.nameOfProfFile,        argumentVector[i+1]);
                }
            i++;
        }
        if (strcmp(argumentVector[i], "-outLog") == 0){
            strcpy(parameter.nameOfLogFile,argumentVector[i+1]);
            i++;
        }
        if (strcmp(argumentVector[i], "-outPress") == 0){
            strcpy(parameter.nameOfDesignationFileForWritingPressureFile,argumentVector[i+1]);
            i++;
        }
        if (strcmp( argumentVector[i], "-outVTK" ) == 0 ) {
            if( parameter.flagOfDivisionOfVtkFile == ON){
                strcpy(parameter.nameOfDividedVtkFile, argumentVector[i+1]);
            }else{
                strcpy(parameter.nameOfVtkFile,        argumentVector[i+1]);
                }
            i++;
        }
  }

/*
  if(argumentCount <= 3) return;

  if( parameter.flagOfDivisionOfProfFile == ON){
	strcpy(parameter.nameOfDividedProfFile, argumentVector[3]);
  }else{
	strcpy(parameter.nameOfProfFile,        argumentVector[3]);
  }

  if(argumentCount <= 4) return;

  strcpy(parameter.nameOfLogFile, argumentVector[4]);

  if(argumentCount <= 5) return;

  strcpy(parameter.nameOfDesignationFileForWritingPressureFile, argumentVector[5]);


  if(argumentCount <= 7) return;

  if( parameter.flagOfDivisionOfVtkFile == ON){
     strcpy(parameter.nameOfDividedVtkFile, argumentVector[7]);
  }else{
     strcpy(parameter.nameOfVtkFile,        argumentVector[7]);
  }

    if(argumentCount <= 8) return;

    //if( parameter.flagOfDivisionOfTorqueFile == ON){
        //strcpy(parameter.nameOfOutputTorqueFile_divided, argumentVector[8]);
    //}else{
        //strcpy(parameter.nameOfOutputTorqueFile,        argumentVector[8]);
    //}
*/
}




void
FILE_displayArgumentOfCommandLine( int argumentCount, char **argumentVector ){

  int iArgument;

  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        Command line arguments                      \n");
  fprintf(FpForLog,"====================================================\n");

  fprintf(FpForLog,"argumentCount = %d\n", argumentCount );

  for( iArgument=0; iArgument < argumentCount; iArgument++){
    fprintf(FpForLog,"argumentVector[%d] = %s\n", iArgument, argumentVector[0] );
  }

}


void
FILE_setDefaultNameOfInputFile( void ){

  strcpy( parameter.nameOfDataFile, NAME_OF_DATA_FILE);

  strcpy( parameter.nameOfGridFile, NAME_OF_GRID_FILE);

  strcpy( parameter.nameOfBubbleInputFile, NAME_OF_BUBBLE_FILE);

  strcpy(parameter.nameOfConcentrationInputFile, NAME_OF_CONCENTRATION_FILE);

}





void
FILE_setDefaultNameOfOutputFile( void ){

  char pathFile[SIZE_OF_FILE_NAME];

  strcpy( parameter.nameOfProfFile, NAME_OF_PROF_FILE);

  strcpy( parameter.nameOfDividedProfFile, NAME_OF_DIVIDED_PROF_FILE);

  strcpy( parameter.nameOfLogFile,  NAME_OF_LOG_FILE);


  strcpy( parameter.nameOfVtkFile, NAME_OF_VTK_FILE);

  strcpy( parameter.nameOfDividedVtkFile, NAME_OF_DIVIDED_VTK_FILE);

  strcpy( parameter.nameOfOutputTorqueFile, NAME_OF_TORQUE_FILE);

  strcpy( parameter.nameOfOutputTorqueFile_divided, NAME_OF_DIVIDED_TORQUE_FILE);


  strcpy( parameter.nameOfMixingIndexFile, NAME_OF_MIXING_INDEX_FILE);


}





void
FILE_readGridFile( void ){

  FILE *fp1;

  FILE *fp2 = NULL;

  FILE *fp3 = NULL;

  int iParticle;
  int iDim;
  int iType;

  printf("\n\nGridFile: %s, ConcFile: %s\n\n",parameter.nameOfGridFile,parameter.nameOfConcentrationInputFile);

  fp1 = FILE_openFile(parameter.nameOfGridFile,"r");

  FILE_scanDouble( fp1,&timer.simulationTime,"timer.simulationTime");
  FILE_scanInt(    fp1,&particle.totalNumber,  "particle.totalNumber");

  particle.totalNumber_upperLimit =  particle.totalNumber;
  fprintf( stderr, "particle.totalNumber_upperLimit = %d\n", particle.totalNumber_upperLimit );

  fprintf(FpForLog,"timer.simulationTime = %lf [sec]\n",       timer.simulationTime);
  fprintf(FpForLog,"particle.totalNumber     = %d  [particles]\n", particle.totalNumber     );

  if(parameter.flagOfBubbleCalculation==ON){
     fp2 = FILE_openFile(parameter.nameOfBubbleInputFile,"r");
  }

  if (parameter.flagOfConcentrationCalculation == ON){
	  fp3 = FILE_openFile(parameter.nameOfConcentrationInputFile, "r");
  }

  MEMORY_allocateMemoryForParticleStructure();

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    FILE_scanInt( fp1,&particle.type[iParticle],"particle.type[iParticle]");

    for(iDim = 0; iDim < 3; iDim++){
      FILE_scanDouble( fp1,&particle.position[iDim][iParticle],"particle.position[iDim][iParticle]");
    }

    for(iDim = 0; iDim < 3; iDim++){
      FILE_scanDouble( fp1,&particle.velocity[iDim][iParticle],"particle.velocity[iDim][iParticle]");
    }

    FILE_scanDouble( fp1,&particle.pressure[iParticle],"particle.pressure[iParticle]");
    FILE_scanDouble( fp1,&particle.particleNumberDensity[iParticle],"particle.particleNumberDensity[iParticle]");

    if(parameter.flagOfBubbleCalculation==ON){
        FILE_scanDouble( fp2,&particle.moleOfBubbles[iParticle],"particle.moleOfBubbles[iParticle]");
        FILE_scanDouble( fp2,&particle.numberOfBubbles[iParticle],"particle.numberOfBubbles[iParticle]");
        FILE_scanDouble( fp2,&particle.concentrationOfImpurities[iParticle],"particle.concentrationOfImpurities[iParticle]");
    }

	if (parameter.flagOfConcentrationCalculation == ON){
		for (iType = 0; iType < 3; iType++){
			if (!strcmp(parameter.nameOfGridFile,parameter.nameOfConcentrationInputFile))
				FILE_scanDouble( fp1, &particle.concentration[iType][iParticle], "particle.concentration[iType][iParticle]");
			else
				FILE_scanDouble( fp3, &particle.concentration[iType][iParticle], "particle.concentration[iType][iParticle]");
		}
	}
  }

  FILE_closeFile( fp1, parameter.nameOfGridFile );

}





void
FILE_readDesignationFileForWritingPressureFile( void ){

  FILE  *fp;
  int   iDesignatedParticle;


  if( parameter.flagOfOutputOfPressureFile == OFF )return;


  fp = FILE_openFile(parameter.nameOfDesignationFileForWritingPressureFile,"r");


  FILE_scanInt(  fp,  &parameter.numberOfDesignatedParticles,    "parameter.numberOfDesignatedParticles");

  fprintf(FpForLog,"parameter.numberOfDesignatedParticles  = %d  [particles]\n", parameter.numberOfDesignatedParticles);


  particle.listOfDesignatedParticles = MEMORY_allocateMemoryFor1dimIntArray(
                                                                       parameter.numberOfDesignatedParticles
                                                                     ,"particle.listOfDesignatedParticles" );


  for(iDesignatedParticle=0; iDesignatedParticle < parameter.numberOfDesignatedParticles; iDesignatedParticle++){

    FILE_scanInt( fp,&particle.listOfDesignatedParticles[iDesignatedParticle]
		      ,"particle.listOfDesignatedParticles[iDesignatedParticle]");
  }


  FILE_closeFile( fp, parameter.nameOfDesignationFileForWritingPressureFile );

}





void
FILE_countNumberOfParticlesEachType( void ){

  int iParticle;

  parameter.numberOfFluidParticles     = 0;
  parameter.numberOfRigidParticles     = 0;
  parameter.numberOfWallParticles      = 0;
  parameter.numberOfDummyWallParticles = 0;
  parameter.numberOfGhostParticles     = 0;
  parameter.numberOfLiquidParticles[0] = 0;
  parameter.numberOfLiquidParticles[1] = 0;
  parameter.numberOfLiquidParticles[2] = 0;
  //printf("wallType: %d, dummyWallType: %d, RigidBodyType: %d, GhostType: %d", parameter.wallType,parameter.dummyWallType,parameter.typeNumberOfRigidParticle_forForcedMotion,GHOST);

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == parameter.wallType ){

      parameter.numberOfWallParticles++;

    }else if( particle.type[iParticle] == parameter.dummyWallType ){

      parameter.numberOfDummyWallParticles++;

    }else if( particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ){

      parameter.numberOfRigidParticles++;

    }else if( particle.type[iParticle] == GHOST ){

      parameter.numberOfGhostParticles++;

    }else{

      parameter.numberOfFluidParticles++;

	  if (particle.type[iParticle] == 0) {
		  parameter.numberOfLiquidParticles[0] += 1;
	  }
	  if (particle.type[iParticle] == 1) parameter.numberOfLiquidParticles[1]++;
	  if (particle.type[iParticle] == 2) parameter.numberOfLiquidParticles[2]++;

    }

  }

}



void
FILE_displayNumberOfParticles( void ){

  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        the number of particles                     \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog, "total                   = %d\n",particle.totalNumber                 );
  fprintf(FpForLog, "  ---fluidParticles     = %d\n",parameter.numberOfFluidParticles     );
  fprintf(FpForLog, "  ---rigidParticles     = %d\n",parameter.numberOfRigidParticles     );
  fprintf(FpForLog, "  ---wallParticles      = %d\n",parameter.numberOfWallParticles      );
  fprintf(FpForLog, "  ---dummyWallParticles = %d\n",parameter.numberOfDummyWallParticles );
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");

}



void
FILE_scanDouble( FILE *fp, double *scanedValue, char *parameterName){

  int status;
  status = fscanf(fp, "%lf", scanedValue);
  FILE_checkEndOfFile(status, parameterName);

}



void
FILE_scanInt( FILE *fp, int *scanedValue, char *parameterName){

  int status;
  status = fscanf(fp, "%d", scanedValue);
  FILE_checkEndOfFile(status, parameterName);

}




void
FILE_scanChar( FILE *fp, char *scanedValue, char *parameterName){

  int status;
  status = fscanf(fp, "%s", scanedValue);
  FILE_checkEndOfFile(status, parameterName);

}



void
FILE_scanOnOff(FILE *fp, int *flag, char *variable_name){
  char buf[256];
  char errorMessage[256];

  fscanf(fp,"%s", buf);
  if( ((strcmp( buf, "ON" ) == 0 )
       ||(strcmp( buf, "on" ) == 0 ))
       ||(strcmp( buf, "On" ) == 0 )){

    (*flag) = ON;

  }else if( ((strcmp( buf, "OFF" ) == 0 )
	     || (strcmp( buf, "off" ) == 0 ))
	     || (strcmp( buf, "Off" ) == 0 )){

    (*flag) = OFF;

  }else{
    sprintf(errorMessage,"Error: Parameter of \"%s\" is not adequate. [ in data-file]\n", variable_name);
    OTHER_endProgram(errorMessage);
  }

}





void
FILE_scanTypeOfOutputInterval(FILE *fp, int *flag, char *variable_name){
  char buf[256];
  char errorMessage[256];

  fscanf(fp,"%s", buf);
  if( ((strcmp( buf, "simulationTime" ) == 0 )
       ||(strcmp( buf, "SimulationTime" ) == 0 ))
       ||(strcmp( buf, "SIMULATION_TIME" ) == 0 )){

    (*flag) = SIMULATION_TIME;


  }else if( ((strcmp( buf, "timeStep" ) == 0 )
	     || (strcmp( buf, "TimeStep" ) == 0 ))
	     || (strcmp( buf, "TIME_STEP" ) == 0 )){


    (*flag) = TIME_STEP;


  }else{
    sprintf(errorMessage,"Error: Parameter of \"%s\" is not adequate. [ in data-file]\n", variable_name);
    OTHER_endProgram(errorMessage);
  }

}



void
FILE_scanMotionTypeOfRigidBody(FILE *fp, int *flag, char *variable_name){
  char buf[256];
  char errorMessage[256];

  fscanf(fp,"%s", buf);
  if( (((strcmp( buf, "forcedMotion" ) == 0 )
	||(strcmp( buf, "forcedmotion" ) == 0 )))
        ||((strcmp( buf, "ForcedMotion" ) == 0 )
        ||(strcmp( buf, "FORCED_MOTION" ) == 0 ))){

    (*flag) = FORCED_MOTION;

  }else if( (((strcmp( buf, "compulsiveMotion" ) == 0 )
	||(strcmp( buf, "compulsivemotion" ) == 0 )))
        ||((strcmp( buf, "CompulsiveMotion" ) == 0 )
        ||(strcmp( buf, "COMPULSIVE_MOTION" ) == 0 ))){

    (*flag) = FORCED_MOTION;

  }else{
    sprintf(errorMessage,"Error: Parameter of \"%s\" is not adequate. [ in data-file]\n", variable_name);
    OTHER_endProgram(errorMessage);
  }

}





char*
FILE_returnOnOff(int flag){

  if(flag == ON){

    return("on");

  }else if(flag == OFF){

    return("off");

  }else{

    return("***--ERROR--****");

  }

}





char*
FILE_returnDim(int iDim){

  if(iDim == XDIM){
    return("XDIM");
  }else if(iDim == YDIM){
    return("YDIM");
  }else if(iDim == ZDIM){
    return("ZDIM");
  }

  return("##---ERROR---##");

}



char*
FILE_returnTypeOfOutputInterval(int flag){

  if(flag == SIMULATION_TIME){

    return("simulationTime");

  }else if(flag == TIME_STEP){

    return("timeStep");

  }else{

    return("***--ERROR--****");

  }

}





char*
FILE_returnMotionTypeOfRigidBody(int flag){

  if(flag == FORCED_MOTION){

    return("forcedMotion");

  }else{

    return("***--ERROR--****");

  }

}



void
FILE_skipComment( FILE *fp ){

  char comment[256];

  fscanf(fp,"%s",comment);

}




void
FILE_checkEndOfFile(int status, char *name){

  if(status == EOF){
    fprintf(FpForLog,"ERROR: Scan of '%s' failed.\n", name);
    fprintf(FpForLog,"\n");
    exit(1);
  }

}


void
FILE_displayGridInformation( void ){

  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        grid information                            \n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"physical time       =  %lf\n", timer.simulationTime);
  fprintf(FpForLog,"number of particles =  %d\n",  particle.totalNumber     );
  fprintf(FpForLog,"\n\n");

}



void
FILE_writeCalculationResultInFile( void ){


  if( YES == TIMER_checkWhetherItIsTimeToWriteProfFile()){

    if((timer.flagOfDisplayingStateOfTimeStep == OFF)&&(parameter.flagOfHideOutput == OFF)){
      TIMER_displayStateOfTimeStep( stderr );
    }

    FILE_writeProfFile();

  }else if(timer.flagOfDisplayingStateOfTimeStep == ON){
    fprintf(stderr," \n");
  }


  FILE_writePressureFile();

}





void
FILE_writeProfFile( void ){

  FILE *fp;
  int   iParticle;
  int   iDim;
  char  fileName[256];
  int   iType;


  fp = FILE_openProfFile();


  fprintf(fp, "%lf\n", timer.simulationTime);
  fprintf(fp, "%d\n",  particle.totalNumber);


  if( parameter.flagOfExponentialExpressionInProfFile == ON ){

	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

	  fprintf(fp,"%d",particle.type[iParticle]);

	  for(iDim = 0; iDim < 3; iDim++){
		fprintf(fp," %e",particle.position[iDim][iParticle]);
	  }

	  for(iDim = 0; iDim < 3; iDim++){
		fprintf(fp," %e",particle.velocity[iDim][iParticle]);
	  }

	  fprintf(fp," %e",particle.pressure[iParticle]);
	  fprintf(fp," %e",particle.particleNumberDensity[iParticle]);

      if( parameter.flagOfBubbleCalculation == ON){
          fprintf(fp," %e",particle.moleOfBubbles[iParticle]);
          fprintf(fp," %e",particle.numberOfBubbles[iParticle]);
          fprintf(fp," %e",particle.concentrationOfImpurities[iParticle]);
      }

	  if (parameter.flagOfConcentrationCalculation == ON){
		  for (iType = 0; iType < 3; iType++){
			  fprintf(fp, " %e", particle.concentration[iType][iParticle]);
		  }
	  }


	  fprintf(fp,"\n");
	}


  }else{

	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
	  fprintf(fp,"%d",particle.type[iParticle]);

	  for(iDim = 0; iDim < 3; iDim++){
		fprintf(fp," %lf",particle.position[iDim][iParticle]);
	  }

	  for(iDim = 0; iDim < 3; iDim++){
		fprintf(fp," %lf",particle.velocity[iDim][iParticle]);
	  }

	  fprintf(fp," %lf",particle.pressure[iParticle]);
	  fprintf(fp," %lf",particle.particleNumberDensity[iParticle]);

      if( parameter.flagOfBubbleCalculation == ON){
            fprintf(fp," %e",particle.moleOfBubbles[iParticle]);
            fprintf(fp," %e",particle.numberOfBubbles[iParticle]);
            fprintf(fp," %e",particle.concentrationOfImpurities[iParticle]);
      }

	  if (parameter.flagOfConcentrationCalculation == ON){
		  for (iType = 0; iType < 3; iType++){
			  fprintf(fp, " %e", particle.concentration[iType][iParticle]);
		  }
	  }

	  fprintf(fp,"\n");
	}

  }



  fflush(fp);
  FILE_closeFile(fp, fileName );

  if(parameter.flagOfHideOutput == OFF) fprintf(stderr,   "%d-th profFile was written.\n", timer.iProfFile);
  fprintf(FpForLog, "%d-th profFile was written.\n", timer.iProfFile);

  FILE_compressProfFile();

  FILE_makeVtkFile();


  timer.iProfFile++;

}


void
FILE_makeVtkFile( void ){

    FILE *fp;
    int   iParticle;
    int   iDim;
	int   iType;
    char filename[256];

    double diameter;
    double Force[3];
    double Advection[3];
    double Torque[3];
    double Y_Torque;

    fp=FILE_openVtkFile();

    //sprintf(filename, "output_%04d.vtk",particle.totalNumber);
    //printf("Creating %s ... ", filename);

    //fp=fopen(filename,"w");
    fprintf(fp,"# vtk DataFile Version 3.0\n");
    fprintf(fp,"vtk output\n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp,"POINTS %d float\n",particle.totalNumber);
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
        for(iDim = 0; iDim < 3; iDim++){
            fprintf(fp,"%lf ",particle.position[iDim][iParticle]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"\n");

    fprintf(fp,"CELLS %d %d\n",particle.totalNumber,particle.totalNumber*2);
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
        fprintf(fp,"1 %d\n",iParticle);
    }
    fprintf(fp,"\n");

    fprintf(fp,"CELL_TYPES %d\n",particle.totalNumber);
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
        fprintf(fp,"1\n");
    }
    fprintf(fp,"\n");

    fprintf(fp,"POINT_DATA %d\n",particle.totalNumber);

    fprintf(fp,"SCALARS type int %d\n",1);
    fprintf(fp,"LOOKUP_TABLE type\n");
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        fprintf(fp,"%d\n",particle.type[iParticle]);
    }
    fprintf(fp,"\n");

    fprintf(fp,"SCALARS Velocity float 1\n");
    fprintf(fp,"LOOKUP_TABLE Velocity\n");
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        double val=sqrt(particle.velocity[0][iParticle]*particle.velocity[0][iParticle]+particle.velocity[1][iParticle]*particle.velocity[1][iParticle]+particle.velocity[2][iParticle]*particle.velocity[2][iParticle]);
        fprintf(fp,"%lf\n",val);
    }
    fprintf(fp,"\n");

    fprintf(fp,"SCALARS Pressure float 1\n");
    fprintf(fp,"LOOKUP_TABLE Pressure\n");
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        fprintf(fp,"%lf\n",particle.pressure[iParticle]);
    }
    fprintf(fp,"\n");

    fprintf(fp,"SCALARS particleNumberDensity float 1\n");
    fprintf(fp,"LOOKUP_TABLE particleNumberDensity\n");
    for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
        fprintf(fp,"%lf\n",particle.particleNumberDensity[iParticle]);
    }
    fprintf(fp,"\n");

	if (parameter.flagOfConcentrationCalculation == ON){
		for (iType = 0; iType < 3; iType++){
			fprintf(fp, "SCALARS ConcentrationLiquid%d float 1\n",iType);
			fprintf(fp, "LOOKUP_TABLE ConcentrationLiquid%d\n", iType);
			for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){
				fprintf(fp, "%lf\n", particle.concentration[iType][iParticle]);
			}
			fprintf(fp, "\n");
		}
	}

	if (parameter.flagOfInterfaceCalculation == ON){
			fprintf(fp, "VECTORS InterfaceGradient float\n");
			for (iParticle = 0; iParticle < particle.totalNumber; iParticle++){
				fprintf(fp,"%lf %lf %lf\n",particle.interfaceGradient[XDIM][iParticle],particle.interfaceGradient[YDIM][iParticle],particle.interfaceGradient[ZDIM][iParticle]);
			}
			fprintf(fp, "\n");
	}



    if(parameter.flagOfBubbleCalculation==ON){

       fprintf(fp,"SCALARS moleOfBubbles double 1\n");
       fprintf(fp,"LOOKUP_TABLE moleOfBubbles\n");
       for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
           fprintf(fp,"%le\n",particle.moleOfBubbles[iParticle]);
       }
       fprintf(fp,"\n");

       fprintf(fp,"SCALARS numberOfBubbles double 1\n");
       fprintf(fp,"LOOKUP_TABLE numberOfBubbles\n");
       for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
           fprintf(fp,"%le\n",particle.numberOfBubbles[iParticle]);
       }
       fprintf(fp,"\n");

       fprintf(fp,"SCALARS RatioOfImpurities double 1\n");
       fprintf(fp,"LOOKUP_TABLE RatioOfImpurities\n");
       for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
           fprintf(fp,"%le\n", particle.concentrationOfImpurities[iParticle]);
       }
       fprintf(fp,"\n");


       fprintf(fp,"SCALARS BubbleDiameter float 1\n");
       fprintf(fp,"LOOKUP_TABLE BubbleDiameter\n");
       for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
           diameter=Bubble_calculateDiameter(iParticle);
           fprintf(fp,"%lf\n", diameter);
       }
       fprintf(fp,"\n");
    }

    fprintf(fp,"VECTORS point_vectors float\n");
	for(iParticle=0;iParticle<particle.totalNumber;iParticle++){
		fprintf(fp,"%lf %lf %lf\n",particle.velocity[XDIM][iParticle],particle.velocity[YDIM][iParticle],particle.velocity[ZDIM][iParticle]);
	}
    fprintf(fp,"\n");

    if(parameter.flagOfOutputOfTorqueFile == ON ){
       fprintf(fp,"VECTORS Torque float\n");

       Y_Torque=0.0;
	   for(iParticle=0;iParticle<particle.totalNumber;iParticle++){

        if(timer.iVtkFile == 0){
            Torque[XDIM]=0.0;
            Torque[YDIM]=0.0;
            Torque[ZDIM]=0.0;
            Y_Torque += 2.0*M_PI*1.0*Torque[YDIM];
            fprintf(fp,"%lf %lf %lf\n",Torque[XDIM],Torque[YDIM],Torque[ZDIM]);
        }else{

           if(particle.type[iParticle]==parameter.typeNumberOfRigidParticle_forForcedMotion){

              for(iDim=0;iDim< NumberOfDimensions;iDim++){
                  Force[iDim]=physicalProperty.massDensity[particle.type[iParticle]]*(4.0/3.0)*M_PI*pow((particle.averageDistance/2.0),3.0)*(1.0/timer.dt)*(particle.velocity3[iDim][iParticle]-particle.velocity2[iDim][iParticle]);
              }
           }else{
              for(iDim=0;iDim< NumberOfDimensions;iDim++){
                   Force[iDim]=0.0;
              }
           }
           Advection[XDIM]=particle.position2[XDIM][iParticle];
           Advection[YDIM]=0.0;
           Advection[ZDIM]=particle.position2[ZDIM][iParticle];

           MATHE_outerProduct(Torque,Advection,Force);
           Y_Torque += 2.0*M_PI*1.0*Torque[YDIM];

           fprintf(fp,"%lf %lf %lf\n",Torque[XDIM],Torque[YDIM],Torque[ZDIM]);
        }
       }
       printf("Torque[Pressure]=%lf\n",Y_Torque);
       fprintf(fp,"\n");
    }

    //printf("done.\n");

    fflush(fp);
    FILE_closeFile(fp, filename);

    timer.iVtkFile++;
}




void
FILE_writePressureFile( void ){

  if( parameter.flagOfOutputOfPressureFile == OFF ) return;

  if( parameter.flagOfOutputOfAllWallParticle == ON ){

    FILE_outputPressureOfAllWallParticles();

  }else{

    FILE_outputPressureOfOnlyDesignatedParticles();

  }

}



void
FILE_outputPressureOfAllWallParticles( void ){

  FILE *fp;
  int  iParticle;
  int  iDim;

  fp = FILE_openPressureFile();

  fprintf(fp, "%lf\n", timer.simulationTime             );
  fprintf(fp, "%d\n",  parameter.numberOfWallParticles );


  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.type[iParticle] == parameter.wallType ){

      fprintf(fp,"%d",iParticle);
      fprintf(fp,"  %d",particle.type[iParticle]);

      for(iDim = 0; iDim < 3; iDim++){
		fprintf(fp," %lf", particle.position[iDim][iParticle]);
      }

      for(iDim = 0; iDim < 3; iDim++){
		fprintf(fp," %lf", particle.velocity[iDim][iParticle]);
      }

      fprintf(fp," %lf", particle.pressure[iParticle]);
      fprintf(fp," %lf", particle.particleNumberDensity[iParticle]);
      fprintf(fp,"\n");

    }

  }


  fflush(fp);
  FILE_closeFile( fp, "output.press" );

  /*
  fprintf(stderr,   "%d-th pressureFile was written.\n", timer.iPressureFile);
  */
  fprintf(FpForLog, "%d-th pressureFile was written.\n", timer.iPressureFile);

  timer.iPressureFile++;

}





void
FILE_outputPressureOfOnlyDesignatedParticles( void ){

  FILE  *fp;
  int   iParticle;
  int   iDesignatedParticle;
  int   iDim;
  char  fileName[256];


  fp = FILE_openPressureFile();

  fprintf(fp, "%lf\n", timer.simulationTime);
  fprintf(fp, "%d\n",  parameter.numberOfDesignatedParticles   );


  for(iDesignatedParticle=0; iDesignatedParticle < parameter.numberOfDesignatedParticles; iDesignatedParticle++){

    iParticle = particle.listOfDesignatedParticles[iDesignatedParticle];

    fprintf(fp,"%d  %d",iParticle, particle.type[iParticle]);

    for(iDim = 0; iDim < 3; iDim++){
      fprintf(fp," %lf",particle.position[iDim][iParticle]);
    }

    for(iDim = 0; iDim < 3; iDim++){
      fprintf(fp," %lf",particle.velocity[iDim][iParticle]);
    }

    fprintf(fp," %lf",particle.pressure[iParticle]);
    fprintf(fp," %lf",particle.particleNumberDensity[iParticle]);

    fprintf(fp,"\n");
  }

  fflush(fp);
  FILE_closeFile(fp, fileName );


  fprintf(FpForLog, "%d-th pressureFile was written.\n", timer.iPressureFile);
  /*
  fprintf(stderr, "%d-th pressureFile was written.\n", timer.iPressureFile);
  */

  timer.iPressureFile++;

}


FILE*
FILE_openProfFile( void ){

  FILE *fp;
  char  fileName[256];


  if( parameter.flagOfDivisionOfProfFile == ON){

     sprintf(fileName,"%s%04d.prof", parameter.nameOfDividedProfFile, timer.iProfFile);
     fp = FILE_openFile(fileName,"w");

  }else{

     sprintf(fileName,"%s", parameter.nameOfProfFile );

     if( timer.iProfFile == timer.initialSequantialNumberOfProfFile ){

       fp = FILE_openFile(fileName,"w");

     }else{

       fp = FILE_openFile(fileName,"a");

     }

  }

  return fp;

}

FILE*
FILE_openVtkFile( void ){

    FILE *fp;
    char  fileName[256];


    if( parameter.flagOfDivisionOfVtkFile == ON){

        sprintf(fileName,"%s%04d.vtk", parameter.nameOfDividedVtkFile, timer.iProfFile);
        fp = FILE_openFile(fileName,"w");

    }else{

        sprintf(fileName,"%s", parameter.nameOfVtkFile );

        if( timer.iProfFile == timer.initialSequantialNumberOfVtkFile ){

            fp = FILE_openFile(fileName,"w");

        }else{

            fp = FILE_openFile(fileName,"a");

        }

    }

    return fp;

}



FILE*
FILE_openPressureFile( void ){

  FILE *fp;
  char  fileName[256];

  if( parameter.flagOfDivisionOfPressureFile == ON){

     sprintf(fileName,"%s%05d.pressure", parameter.nameOfOutputPressureFile_divided, timer.iPressureFile);
     fp = FILE_openFile(fileName,"w");

  }else{

     sprintf(fileName,"%s", parameter.nameOfOutputPressureFile );

     if( timer.iProfFile == 0 ){

       fp = FILE_openFile(fileName,"w");

     }else{

       fp = FILE_openFile(fileName,"a");

     }
  }

  return fp;

}



void
FILE_compressProfFile( void ){


  char  fileName[SIZE_OF_FILE_NAME];
  char  command[SIZE_OF_FILE_NAME];

  if( parameter.flagOfCompressionOfProfFile != ON) return;

  if( parameter.flagOfDivisionOfProfFile == ON){

     sprintf(fileName,"%s%04d.prof", parameter.nameOfDividedProfFile, timer.iProfFile);
     sprintf(command,"gzip -9 %s",fileName);
     system(command);

     fprintf(stderr,   "%d-th profFile was compressed.\n", timer.iProfFile);
     fprintf(FpForLog, "%d-th profFile was compressed.\n", timer.iProfFile);

  }

}



void
FILE_readDataFile(void){

	int  iType;
	int  iDim;
	char variableName[256];
	FILE *fp;


	fp = FILE_openFile(parameter.nameOfDataFile, "r");

	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#######--FREQUENTLY-USED-DATA-#####################################
	//#--------ParticleSize----------------------------------------------

	FILE_scanDouble(fp, &particle.averageDistance           //scan averageDistanceBetweenParticles
		, "particle.averageDistance");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------Time------------------------------------------------------

	FILE_scanDouble(fp, &timer.finishTime                  //scan finishTime
		, "timer.finishTime");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &timer.dt_initial                  //scan initial Dt
		, "timer.dt_initial");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------File------------------------------------------------------

	FILE_scanTypeOfOutputInterval(fp, &timer.typeOfOutputInterval    //scan Type of output Interval
		, "typeOfOutputInterval");


	FILE_skipComment(fp);
	FILE_scanDouble(fp, &timer.intervalTimeOfWritingProfFile_simulationTime    //scan simulation time
		, "timer.intervalTimeOfWritingProfFile_simulationTime");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &timer.intervalTimeOfWritingProfFile_timeStep     //scan timestep
		, "timer.intervalTimeOfWritingProfFile_timeStep");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------Domain----------------------------------------------------

	FILE_scanInt(fp, &NumberOfDimensions                    //scan Number of Dimensions
		, "NumberOfDimensions");

	FILE_skipComment(fp);


	FILE_scanOnOff(fp, &domain.flagOfAutoSettingOfDomainSize      //scan auto setting of domain size
		, "domain.flagOfAutoSettingOfDomainSize");


	for (iDim = 0; iDim < 3; iDim++){                       //scan upper margin ratio
		FILE_skipComment(fp);
		sprintf(variableName, "domain.upperMarginRatio[%d]", iDim);
		FILE_scanDouble(fp, &domain.upperMarginRatio[iDim], variableName);
	}

	for (iDim = 0; iDim < 3; iDim++){                       //scan lower margin ratio
		FILE_skipComment(fp);
		sprintf(variableName, "domain.lowerMarginRatio[%d]", iDim);
		FILE_scanDouble(fp, &domain.lowerMarginRatio[iDim], variableName);
	}


	for (iDim = 0; iDim < 3; iDim++){                      //scan upper limit
		FILE_skipComment(fp);
		sprintf(variableName, "domain.upperLimit[%d]", iDim);
		FILE_scanDouble(fp, &domain.upperLimit[iDim], variableName);
	}

	for (iDim = 0; iDim < 3; iDim++){                      //scan lower limit
		FILE_skipComment(fp);
		sprintf(variableName, "domain.lowerLimit[%d]", iDim);
		FILE_scanDouble(fp, &domain.lowerLimit[iDim], variableName);
	}


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------RadiusOfInteractionDomain--------------------------------

	FILE_scanDouble(fp, &parameter.radiusOfParticleNumberDensity_ratio      //scan radius of particle number density
		, "parameter.radiusOfParticleNumberDensity_ratio");


	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.radiusOfGradient_ratio                  //scan radius of gradient
		, "parameter.radiusOfGradient_ratio");


	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.radiusOfLaplacianForViscosity_ratio      //scan radius of laplacian for viscosity
		, "parameter.radiusOfLaplacianForViscosity_ratio");


	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.radiusOfLaplacianForPressure_ratio      //scan radius of Laplacian for Pressure
		, "parameter.radiusOfLaplacianForPressure_ratio");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.radiusOfLaplacianForDiffusion_ratio     //scan radius of Laplacian for Pressure
		, "parameter.radiusOfLaplacianForDiffusion_ratio");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#######--LESS-FREQUENTLY-USED-DATA-###############################
	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------Bucket(subDomain)----------------------------------------
	FILE_scanOnOff(fp, &domain.flagOfAutoSettingOfBucketCapacity     //scan autosetting of Bucket Capacity
		, "domain.flagOfAutoSettingOfBucketCapacity");

	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &domain.flagOfOptimizationOfBucketMemorySize     //scan optimization of memorzsize
		, "domain.flagOfOptimizationOfBucketMemorySize");


	FILE_skipComment(fp);
	FILE_scanDouble(fp, &domain.marginRatioForSettingBucketCapacity     //scan marginratio
		, "domain.marginRatioForSettingBucketCapacity");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &domain.capacityOfBucket                   //scan capacity of bucket
		, "domain.capacityOfBucket");


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------TableOfNeighborParticles---------------------------------

	FILE_scanInt(fp, &parameter.numberOfNeighborTables       //scan Number of Neighbortables
		, "parameter.numberOfNeighborTables");


	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfAutoSettingOfCapacityOfNeighborTable     //scan autosetting capacity of neighbortable
		, "parameter.flagOfAutoSettingOfCapacityOfNeighborTable");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.marginRatioForSettingCapacityOfLargeNeighborTable     //scan marginratio of large table
		, "parameter.marginRatioForSettingCapacityOfLargeNeighborTable");


	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.marginRatioForSettingCapacityOfSmallNeighborTable     //scan marginratio of small table
		, "parameter.marginRatioForSettingCapacityOfSmallNeighborTable");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.capacityOfNeighborTable_large          //scan capacity of large neighbortable
		, "parameter.capacityOfNeighborTable_large");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.capacityOfNeighborTable_small         //scan capacity of small neighbortable
		, "parameter.capacityOfNeighborTable_small");


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------ProfFile-------------------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfDivisionOfProfFile           //scan division of prof file
		, "parameter.flagOfDivisionOfProfFile");

	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfCompressionOfProfFile         //scan compression of prof file
		, "parameter.flagOfCompressionOfProfFile");


	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.upperLimitOfNumberOfProfFiles         //scan upper limit of Number of files
		, "parameter.upperLimitOfNumberOfProfFiles");


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------PressureFile---------------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfOutputOfPressureFile        //scan output pressurefile
		, "parameter.flagOfOutputOfPressureFile");

	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfOutputOfAllWallParticle     //scan outputpressurefile of all wall particle
		, "parameter.flagOfOutputOfAllWallParticle");

	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfDesignationFileForWritingPressureFile     //scan name of inputfile
		, "parameter.nameOfDesignationFileForWritingPressureFile");

	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfDivisionOfPressureFile       //scan division of pressurefile
		, "parameter.flagOfDivisionOfPressureFile");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.upperLimitOfNumberOfPressureFiles     //scan upper limit of number of files
		, "parameter.upperLimitOfNumberOfPressureFiles");


	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfOutputPressureFile_divided     //scan name of outputfile on
		, "parameter.nameOfOutputPressureFile_divided");

	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfOutputPressureFile        //scan name of outputfile off
		, "parameter.nameOfOutputPressureFile");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------TorqueFile---------------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfOutputOfTorqueFile     //scan output Torquefile
		, "parameter.flagOfOutputOfTorqueFile");

	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfOutputTorqueFile         //scan name of Outputfile
		, "parameter.nameOfOutputTorqueFile");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.numberOfRotations            //scan number of rotations
		, "parameter.numberOfRotations");


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------TypeOfParticle-------------------------------------------

	FILE_scanInt(fp, &parameter.numberOfParticleTypes        //scan number of particle types
		, "parameter.numberOfParticleTypes");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.wallType                    //scan number of wall particles
		, "parameter.wallType");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.dummyWallType              //scan number of dummywall particle
		, "parameter.dummyWallType");


	physicalProperty.massDensity = MEMORY_allocateMemoryFor1dimDoubleArray(parameter.numberOfParticleTypes, "massDenisity");

	physicalProperty.compressibility = MEMORY_allocateMemoryFor1dimDoubleArray(parameter.numberOfParticleTypes, "compressibility");


	FILE_skipComment(fp);

	//#--------MassDensity----------------------------------------------

	for (iType = 0; iType < parameter.numberOfParticleTypes; iType++){      //scan mass density
		FILE_skipComment(fp);
		FILE_scanDouble(fp, &physicalProperty.massDensity[iType]
			, "physicalProperty.massDensity[iType]");
	}


	FILE_skipComment(fp);

	//#--------Compressibility------------------------------------------

	for (iType = 0; iType < parameter.numberOfParticleTypes; iType++){     //scan compressibility
		FILE_skipComment(fp);
		FILE_scanDouble(fp, &physicalProperty.compressibility[iType]
			, "physicalProperty.compressibility[iType]");
	}

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------Viscosity------------------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfViscosityCalculation        //scan flag of viscositycalculation
		, "parameter.flagOfViscosityCalculation");

	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfHighViscosityCalculation     //scan flag of high viscositycalculation
		, "parameter.flagOfHighViscosityCalculation");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &physicalProperty.kinematicViscosity        //scan kinematic viscosity
		, "physicalProperty.kinematicViscosity");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------Bubble---------------------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfBubbleCalculation       //scan flag of bubblecalculation
		, "parameter.flagOfBubbleCalculation");

	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfBubbleInputFile         //scan name of inputfile
		, "parameter.nameOfBubbleInputFile");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.numperOfParticleForCalculatingBeta      //scan number of particles for calculating beta
		, "parameter.numperOfParticleForCalculatingBeta");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &physicalProperty.massDensityOfBubble      //scan mass density of bubble
		, "physicalProperty.massDensityOfBubble");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &physicalProperty.gasConstant          //scan gas constant
		, "physicalProperty.gasConstant");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &physicalProperty.temperature        //scan temperature
		, "physicalProperty.temperature");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &physicalProperty.headPressure        //scan headpressure
		, "physicalProperty.headPressure");

	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);

	//#--------Gravity--------------------------------------------------

	for (iDim = 0; iDim < 3; iDim++){                     //scan gravity
		FILE_skipComment(fp);
		sprintf(variableName, "physicalProperty.gravity[%d]", iDim);
		FILE_scanDouble(fp, &physicalProperty.gravity[iDim], variableName);
	}


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------TimeDifference-------------------------------------------

	FILE_scanDouble(fp, &parameter.courantNumber          //scan courant Number
		, "parameter.courantNumber");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.diffusionNumber       //scan diffusionnumber
		, "parameter.diffusionNumber");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &timer.maxDt_ratio               //scan maxDt
		, "timer.maxDt_ratio");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &timer.minDt_ratio              //scan minDt
		, "timer.minDt_ratio");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &timer.upperLimitOfChangeRateOfDt         //scan upper limit of changerate of dt
		, "timer.upperLimitOfChangeRateOfDt");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------SolvelerOfSimultaneousEquations--------------------------

	FILE_scanInt(fp, &parameter.maxIterationNumberInIterationSolver         //scan upper limit of iterationnumber
		, "parameter.maxIterationNumberInIterationSolver");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.minIterationNumberInIterationSolver       //scan lower limit of iterationnumber
		, "parameter.minIterationNumberInIterationSolver");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.smallNumberForCheckingConvergenceInIterationSolver      //scan small number for checking convergence
		, "parameter.smallNumberForCheckingConvergenceInIterationSolver");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------CollisionBetweenParticles--------------------------------

	FILE_scanDouble(fp, &parameter.collisionDistance_ratio           //scan collision distance
		, "parameter.collisionDistance_ratio");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.collisionCoefficient             //scan collision coefficient
		, "parameter.collisionCoefficient");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------DirichletBoundaryCondition-------------------------------

	FILE_scanDouble(fp, &parameter.thresholdOfParticleNumberDensity_ratio     //scan threshold of particle number density ratio
		, "parameter.thresholdOfParticleNumberDensity_ratio");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------Other----------------------------------------------------

	FILE_scanDouble(fp, &parameter.relaxationCoefficientOfLambda      //scan relaxation coefficient of lambda
		, "parameter.relaxationCoefficientOfLambda");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &timer.finishTimeStep                  //scan finish time step
		, "timer.finishTimeStep");

	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//######---CustomizedFunctions---###################################
	//#--------RigidBody------------------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfForcedMotionOfRigidBody       //scan forced motion of rigid body
		, "parameter.flagOfForcedMotionOfRigidBody");

	FILE_skipComment(fp);
	FILE_scanInt(fp, &parameter.typeNumberOfRigidParticle_forForcedMotion      //scan type number of rigid particle
		, "parameter.typeNumberOfRigidParticle_forForcedMotion");

	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody      //scan filename of samplingdata for rigid body
		, "parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------TANAKAandMASUNAGAModel-----------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfTanakaAndMasunagaModel       //scan flag of tanaka and masunaga model
		, "parameter.flagOfTanakaAndMasunagaModel");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.valueOfGamma                   //scan value of gamma
		, "parameter.valueOfGamma");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.valueOfC                        //scan value of C
		, "parameter.valueOfC");

	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------ConcentrationCalculation---------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfConcentrationCalculation       //scan flag of Concentration Calculation
		, "parameter.flagOfConcentrationCalculation");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.DiffusionCoefficient                 //scan Diffusion Coefficient
		, "parameter.DiffusionCoefficient");

	FILE_skipComment(fp);
	FILE_scanChar(fp, parameter.nameOfConcentrationInputFile      //scan filename of inputfile
		, "parameter.nameOfConcentrationInputFile ");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------InterfaceCalculation--------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfInterfaceCalculation       //scan flag of Interface Calculation
		, "parameter.flagOfInterfaceCalculation");

	FILE_scanChar(fp, variableName, "InterfaceIntervalOrCorrection");
	if(!strcmp(variableName,"intervalTimeOfInterface")){
		FILE_scanDouble(fp, &timer.intervalTimeOfInterface    //scan Interface Interval time
			, "timer.intervalTimeOfInterface");
		FILE_skipComment(fp);
	}
	else
		timer.intervalTimeOfInterface = timer.dt_initial;
	
	/*
	FILE_skipComment(fp);
	FILE_scanDouble(fp, &timer.intervalTimeOfInterface    //scan Interface Interval time
		, "timer.intervalTimeOfInterface");
	
	FILE_skipComment(fp);*/
	FILE_scanOnOff(fp, &parameter.flagOfInterfaceCorrection       //scan flag of Interface Correction
		, "parameter.flagOfInterfaceCorrection");
	
	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfInterfaceCompare          //scan flag of Interface Compare
		, "parameter.flagOfInterfaceCompare");

	FILE_skipComment(fp);
	FILE_scanOnOff(fp, &parameter.flagOfPeriodicX       //scan flag of Periodic Boundary in X direction
		, "parameter.flagOfPeriodicX");

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	//#--------ParticleRedistribution------------------------------------

	FILE_scanOnOff(fp, &parameter.flagOfSingleVortex       //scan flag of Single Vortex Calculation
		, "parameter.flagOfSingleVortex");

	FILE_skipComment(fp);
	FILE_scanDouble(fp, &parameter.alphaZero       //scan flag of Particle Redistribution coefficient alpha
		, "parameter.alphaZero");

	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);
	FILE_skipComment(fp);



	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);


	FILE_skipComment(fp);
	FILE_skipComment(fp);

	FILE_skipComment(fp);
	FILE_skipComment(fp);


	parameter.increaseRateOfDiagonalTermInCoefficientMatrixOfPressure = 2.0;

	parameter.flagOfExponentialExpressionInProfFile = OFF;


	parameter.flagOfDivisionOfVtkFile = ON;


	timer.intervalTimeOfDisplayingStateOfTimeStep = 0.5;


	FILE_closeFile(fp, parameter.nameOfDataFile);

}


void
FILE_displayReadDataFile( void ){

  int  iType;
  int  iDim;

  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"\n");
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"        %s                                 \n",parameter.nameOfDataFile);
  fprintf(FpForLog,"====================================================\n");
  fprintf(FpForLog,"\n");

  fprintf(FpForLog,"#######--FREQUENTLY-USED-DATA-#####################################\n");
  fprintf(FpForLog,"#--------ParticleSize----------------------------------------------\n");
  fprintf(FpForLog,"averageDistanceBetweenParticles(m)                 %lf\n",particle.averageDistance);
  fprintf(FpForLog,"#--------Time------------------------------------------------------\n");
  fprintf(FpForLog,"finishTime(sec)                                    %lf\n",timer.finishTime);
  fprintf(FpForLog,"inititalDt(sec)                                    %lf\n",timer.dt_initial);
  fprintf(FpForLog,"#--------File------------------------------------------------------\n");
  fprintf(FpForLog,"TypeOfOutputInterval(simulationTime/timeStep)      %s\n",FILE_returnTypeOfOutputInterval(timer.typeOfOutputInterval));
  fprintf(FpForLog,"--simulationTime---interval(sec)                   %lf\n",timer.intervalTimeOfWritingProfFile_simulationTime);
  fprintf(FpForLog,"--timeStep---------interval(timeStep)              %d\n",timer.intervalTimeOfWritingProfFile_timeStep);
  fprintf(FpForLog,"#--------Domain----------------------------------------------------\n");
  fprintf(FpForLog,"numberOfDimensions(2or3)                           %d\n",NumberOfDimensions);
  fprintf(FpForLog,"autoSettingOfDomainSize(on/off)                    %s\n",FILE_returnOnOff(domain.flagOfAutoSettingOfDomainSize));

  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"--on---upperMarginRatio[%s](ratio)               %lf\n",FILE_returnDim(iDim),domain.upperMarginRatio[iDim]);
  }

  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"--on---lowerMarginRatio[%s](ratio)               %lf\n",FILE_returnDim(iDim),domain.lowerMarginRatio[iDim]);
  }

  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"--off--upperLimit[%s](m)                         %lf\n",FILE_returnDim(iDim),domain.upperLimit[iDim]);
  }

  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"--off--lowerLimit[%s](m)                         %lf\n",FILE_returnDim(iDim),domain.lowerLimit[iDim]);
  }

  fprintf(FpForLog,"#--------RadiusOfInteractionDomain--------------------------------\n");
  fprintf(FpForLog,"radiusOfParticleNumberDensity(ratio)               %lf\n",parameter.radiusOfParticleNumberDensity_ratio);
  fprintf(FpForLog,"radiusOfGradient(ratio)                            %lf\n",parameter.radiusOfGradient_ratio);
  fprintf(FpForLog,"radiusOfLaplacianForViscosity(ratio)               %lf\n",parameter.radiusOfLaplacianForViscosity_ratio);
  fprintf(FpForLog,"radiusOfLaplacianForPressure(ratio)                %lf\n",parameter.radiusOfLaplacianForPressure_ratio);
  fprintf(FpForLog, "radiusOfLaplacianForDiffusion(ratio)              %lf\n", parameter.radiusOfLaplacianForDiffusion_ratio);
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#######--LESS-FREQUENTLY-USED-DATA-###############################\n");
  fprintf(FpForLog,"#--------Bucket(subDomain)----------------------------------------\n");
  fprintf(FpForLog,"autoSettingOfBucketCapacity                        %s\n",FILE_returnOnOff(domain.flagOfAutoSettingOfBucketCapacity));
  fprintf(FpForLog,"--on---optimizationOfMemorySize(on/off)            %s\n",FILE_returnOnOff(domain.flagOfOptimizationOfBucketMemorySize));
  fprintf(FpForLog,"--on-----off--marginRatio(ratio)                   %lf\n",domain.marginRatioForSettingBucketCapacity);
  fprintf(FpForLog,"--off--capacityOfBucket                            %d\n",domain.capacityOfBucket);
  fprintf(FpForLog,"#--------TableOfNeighborParticles---------------------------------\n");
  fprintf(FpForLog,"theNumberOfNeighborTables(1or2)                    %d\n",parameter.numberOfNeighborTables);
  fprintf(FpForLog,"autoSettingOfCapacityOfNeighborTable               %s\n",FILE_returnOnOff(parameter.flagOfAutoSettingOfCapacityOfNeighborTable));
  fprintf(FpForLog,"--on---marginRatioOfLargeTable(ratio)              %lf\n",parameter.marginRatioForSettingCapacityOfLargeNeighborTable);
  fprintf(FpForLog,"--on---marginRatioOfSmallTable(ratio)              %lf\n",parameter.marginRatioForSettingCapacityOfSmallNeighborTable);
  fprintf(FpForLog,"--off--capacityOfLargeNeighborTable                %d\n",parameter.capacityOfNeighborTable_large);
  fprintf(FpForLog,"--off--capacityOfSmallNeighborTable                %d\n",parameter.capacityOfNeighborTable_small);

  fprintf(FpForLog,"#--------ProfFile-------------------------------------------------\n");
  fprintf(FpForLog,"divisionOfProfFile(on/off)                         %s\n",FILE_returnOnOff(parameter.flagOfDivisionOfProfFile));
  fprintf(FpForLog,"--on---compressionOfProfFile                       %s\n",FILE_returnOnOff(parameter.flagOfCompressionOfProfFile));
  fprintf(FpForLog,"--on---upperLimitOfNumberOfFiles                   %d\n",parameter.upperLimitOfNumberOfProfFiles);

  fprintf(FpForLog,"#--------PressureFile---------------------------------------------\n");
  fprintf(FpForLog,"outputPressureFile(on/off)                         %s\n",FILE_returnOnOff(parameter.flagOfOutputOfPressureFile));
  fprintf(FpForLog,"--on---outputPressureOfAllWallParticle(on/off)     %s\n",FILE_returnOnOff(parameter.flagOfOutputOfAllWallParticle));
  fprintf(FpForLog,"---------off--nameOfInputFile                      %s\n",parameter.nameOfDesignationFileForWritingPressureFile);
  fprintf(FpForLog,"--on---divisionOfPressureFile(on/off)              %s\n",FILE_returnOnOff(parameter.flagOfDivisionOfPressureFile));
  fprintf(FpForLog,"---------on---upperLimitOfNumberOfFiles            %d\n",parameter.upperLimitOfNumberOfPressureFiles);
  fprintf(FpForLog,"---------on---nameOfOutputFile                     %s\n",parameter.nameOfOutputPressureFile_divided);
  fprintf(FpForLog,"---------on---nameOfOutputFile                     %s\n",parameter.nameOfOutputPressureFile);

  fprintf(FpForLog,"#--------TorqueFile---------------------------------------------\n");
  fprintf(FpForLog,"outputTorqueFile(on/off)                         %s\n",FILE_returnOnOff(parameter.flagOfOutputOfTorqueFile));
  //fprintf(FpForLog,"--on---outputTorqueOfAllRigidParticle(on/off)     %s\n",FILE_returnOnOff(parameter.flagOfOutputOfAllRigidParticle));
  //fprintf(FpForLog,"---------off--nameOfInputFile                      %s\n",parameter.nameOfDesignationFileForWritingTorqueFile);
  //fprintf(FpForLog,"--on---divisionOfTorqueFile(on/off)              %s\n",FILE_returnOnOff(parameter.flagOfDivisionOfTorqueFile));
  //fprintf(FpForLog,"---------on---upperLimitOfNumberOfFiles            %d\n",parameter.upperLimitOfNumberOfTorqueFiles);
  //fprintf(FpForLog,"---------on---nameOfOutputFile                     %s\n",parameter.nameOfOutputTorqueFile_divided);
  fprintf(FpForLog,"---------on---nameOfOutputFile                     %s\n",parameter.nameOfOutputTorqueFile);

  fprintf(FpForLog,"#--------TypeOfParticle-------------------------------------------\n");
  fprintf(FpForLog,"numberOfParticleTypes                              %d\n",parameter.numberOfParticleTypes);
  fprintf(FpForLog,"--typeNumberOfWallParticle                         %d\n",parameter.wallType);
  fprintf(FpForLog,"--typeNumberOfDummyWallParticle                    %d\n",parameter.dummyWallType);

  fprintf(FpForLog,"#--------MassDensity----------------------------------------------\n");
  for(iType=0; iType < parameter.numberOfParticleTypes; iType++){
    fprintf(FpForLog,"massDensityOfType%d(kg/m^3)                         %lf\n",iType, physicalProperty.massDensity[iType]);
  }

  fprintf(FpForLog,"#--------Compressibility------------------------------------------\n");
  for(iType=0; iType < parameter.numberOfParticleTypes; iType++){
    fprintf(FpForLog,"compressibilityOfType%d                             %e\n",iType, physicalProperty.compressibility[iType]);
  }

  fprintf(FpForLog,"#--------Viscosity------------------------------------------------\n");
  fprintf(FpForLog,"flagOfViscosityCalculation(on/off)                 %s\n", FILE_returnOnOff(parameter.flagOfViscosityCalculation));

  fprintf(FpForLog,"--on--flagOfHighViscosityCalculation(on/off)       %s\n", FILE_returnOnOff(parameter.flagOfHighViscosityCalculation));

  fprintf(FpForLog,"--on--kinematicViscosity(m^2/s)                    %e\n", physicalProperty.kinematicViscosity);

  fprintf(FpForLog,"#--------Bubble------------------------------------------------\n");
  fprintf(FpForLog,"flagOfBubbleCalculationy(on/off)                  %s\n", FILE_returnOnOff(parameter.flagOfBubbleCalculation));
  fprintf(FpForLog,"--on---nameOfInputfile                            %s\n", parameter.nameOfBubbleInputFile);
  fprintf(FpForLog,"--on--massDensityOfBubble(kg/m^3)                 %e\n", physicalProperty.massDensityOfBubble);
  fprintf(FpForLog,"--on--gassConstant(J/K*mol)                       %e\n", physicalProperty.temperature);
  fprintf(FpForLog,"--on--temperature(K)                              %e\n", physicalProperty.gasConstant);
  fprintf(FpForLog,"--on--headPressure(Pa)                            %e\n", physicalProperty.headPressure);

  fprintf(FpForLog,"#-----------------------------------------------------------------\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");

  fprintf(FpForLog,"#--------Gravity--------------------------------------------------\n");
  for(iDim=0; iDim < 3; iDim++){
    fprintf(FpForLog,"gravity[%s](m/s^2)                               %lf\n",FILE_returnDim(iDim),physicalProperty.gravity[iDim]);
  }

  fprintf(FpForLog,"#--------TimeDifference-------------------------------------------\n");
  fprintf(FpForLog,"courantNumber                                      %lf\n",parameter.courantNumber);
  fprintf(FpForLog,"diffusionNumber                                    %lf\n",parameter.diffusionNumber);
  fprintf(FpForLog,"maxDt(ratio)                                       %lf\n",timer.maxDt_ratio);
  fprintf(FpForLog,"minDt(ratio)                                       %lf\n",timer.minDt_ratio);
  fprintf(FpForLog,"upperLimitOfChangeRateOfDt(ratio)                  %lf\n",timer.upperLimitOfChangeRateOfDt);

  fprintf(FpForLog,"#--------SolvelerOfSimultaneousEquations--------------------------\n");
  fprintf(FpForLog,"upperLimitOfIterationNumber                        %d\n",parameter.maxIterationNumberInIterationSolver);
  fprintf(FpForLog,"lowerLimitOfIterationNumber                        %d\n",parameter.minIterationNumberInIterationSolver);
  fprintf(FpForLog,"smallNumberForCheckingConvergence                  %e\n",parameter.smallNumberForCheckingConvergenceInIterationSolver);

  fprintf(FpForLog,"#--------CollisionBetweenParticles--------------------------------\n");
  fprintf(FpForLog,"collisionDistance(ratio)                           %lf\n",parameter.collisionDistance_ratio);
  fprintf(FpForLog,"collisionCoefficient                               %lf\n",parameter.collisionCoefficient);

  fprintf(FpForLog,"#--------DirichletBoundaryCondition-------------------------------\n");
  fprintf(FpForLog,"thresholdOfParticleNumberDensity(ratio)            %lf\n",parameter.thresholdOfParticleNumberDensity_ratio);


  fprintf(FpForLog,"#--------Other----------------------------------------------------\n");
  fprintf(FpForLog,"relaxationCoefficientOfLambda                        %lf\n",parameter.relaxationCoefficientOfLambda);
  fprintf(FpForLog,"finishTimeStep(timeStep)                           %d\n",timer.finishTimeStep);


  fprintf(FpForLog,"######---CustomizedFunctions---###################################\n");
  fprintf(FpForLog,"#--------RigidBody------------------------------------------------\n");
  fprintf(FpForLog,"forcedMotionOfRigidBody(on/off)                    %s\n",FILE_returnOnOff(parameter.flagOfForcedMotionOfRigidBody));
  fprintf(FpForLog,"--on---typeNumberOfRigidParticle                   %d\n",parameter.typeNumberOfRigidParticle_forForcedMotion);
  fprintf(FpForLog,"--on---fileNameOfSamplingDataForRigidBody          %s\n",parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody);
  fprintf(FpForLog,"#--------TANAKAandMASUNAGA------------------------------------------------\n");
    fprintf(FpForLog,"flagOfTanakaAndMasunagaModel(on/off)                    %s\n",FILE_returnOnOff(parameter.flagOfTanakaAndMasunagaModel));
    fprintf(FpForLog,"--on--valueOfGamma                   %lf\n",parameter.valueOfGamma);
    fprintf(FpForLog,"--on--valueOfC          %lf\n",parameter.valueOfC);
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#--------ConcentrationCalculation---------------------------------\n");
  fprintf(FpForLog,"flagOfConcentrationCalculation(on/off)             %s\n", FILE_returnOnOff(parameter.flagOfConcentrationCalculation));
  fprintf(FpForLog,"DiffusionCoefficient                               %lf\n", parameter.DiffusionCoefficient);
  fprintf(FpForLog,"--on---nameofInputfile                             %s\n", parameter.nameOfConcentrationInputFile);
  fprintf(FpForLog,"#-----------------------------------------------------------------\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");

  fprintf(FpForLog,"#-----------------------------------------------------------------\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");

  fprintf(FpForLog,"#-----------------------------------------------------------------\n");
  fprintf(FpForLog,"#############                                      ***\n");

  fprintf(FpForLog,"#-----------------------------------------------------------------\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"#############                                      ***\n");
  fprintf(FpForLog,"\n\n\n");

}

