#include <time.h>



typedef struct{

  int    *type;
  int    *type_initial;

  double  **position;
  double  **position_previous;

  double  **velocity;
  double  **velocity_previous;
  double  **velocity_correction;

  double  *velocity_component;



  double *pressure;
  double *pressure_previous;


  double *particleNumberDensity;
  double *particleNumberDensity_previous;
  double *particleNumberDensity_weightIsOne;
  double *particleNumberDensity_withConstantWeight;

  double **particleNumberDensity_eachParticleType;

  int    *flagOfBoundaryCondition;
  double *sourceTermOfPressure;
  /*
  double *sourceTermOfPressure_forBiCG;
  */

  double *minPressureAmongNeighbors;

  int    totalNumber;

  int    totalNumber_upperLimit;
  int    standardParticleNumber;

  double minPosition[3];
  double minPosition_initial[3];
  double minVelocity[3];
  double minPressure;
  double minParticleNumberDensity;

  double maxPosition[3];
  double maxPosition_initial[3];
  double maxVelocity[3];
  double maxPressure;
  double maxParticleNumberDensity;
  double maxAbsoluteValueOfVelocity;

  double **coefficientMatrixOfPressure;

  double **coefficientMatrixOfViscosity;
  double *sourceTermOfViscosity;

  /*
  double **coefficientMatrixOfPressure_forBiCG;
  */

  int     sizeOfMatrixArray;

  int    **neighborTable_large;
  int    **neighborTable_small;

  int    **neighborTable_largeBoundary;
  int    **neighborTable_smallBoundary;

  int    *numberOfNeighborParticles_large;
  int    *numberOfNeighborParticles_small;

  int    *listOfDesignatedParticles;

  double averageDistance;

  double *moleOfBubbles;
  double *numberOfBubbles;

  double *moleOfBubbles_previous;
  double *numberOfBubbles_previous;

  double *concentrationOfImpurities;
  double *concentrationOfImpurities_previous;

  double **position1;
  double **position2;
  double **position3;
  double **position4;
  double **velocity1;
  double **velocity2;
  double **velocity3;
  double **velocity4;

  double **concentration;
  double **concentration_previous;

  double **interfaceGradient;
  double *interfaceArea;

}structParticle;



typedef struct{

  double dt;
  double dt_initial;

  double minDt;
  double maxDt;

  double minDt_ratio;
  double maxDt_ratio;

  double simulationTime;
  double finishTime;

  int    iTimeStep;
  int    iTimeStep_copy;

  int    finishTimeStep;

  double upperLimitOfChangeRateOfDt;

  int    typeOfOutputInterval;
  double intervalTimeOfWritingProfFile_simulationTime;
  int    intervalTimeOfWritingProfFile_timeStep;

  double intervalTimeOfDisplayingStateOfTimeStep;

  double timeToWriteProfFile_simulationTime;
  int    timeToWriteProfFile_timeStep;
  double timeToDisplayStateOfTimeStep;

  int    iProfFile;
  int    initialSequantialNumberOfProfFile;

  int    initialSequantialNumberOfVtkFile;

  int    iPressureFile;


  //int    iTorqueFile;
  int    iVtkFile;


  clock_t clockAtStartingTime;
  clock_t clockAtCurrentTime;
  clock_t clockAtBeginningOfMainLoop;


  double  calculationTime;
  double  calculationTimePerTimeStep;

  int     hours;
  int     minutes;
  double  seconds;

  int     flagOfDisplayingStateOfTimeStep;

  double realTimeStart;
  double realTimeCurrent;
  double realTimeLoop;

  double realTime;
  double realTimePerTimeStep;
  int    realHours;
  int    realMinutes;
  double realSeconds;

  double intervalTimeOfInterface;

}structTimer;



typedef struct{

  double  gravity[3];
  double *massDensity;
  double *compressibility;
  double  kinematicViscosity;

  double  massDensityOfBubble;
  double  gasConstant;
  double  temperature;
  double  headPressure;

}structPhysicalProperty;



typedef struct{

  double radiusOfParticleNumberDensity;
  double radiusOfParticleNumberDensity_ratio;

  double radiusOfGradient;
  double radiusOfGradient_ratio;
  double radiusOfGradient_squared;
  double radiusOfLaplacianForInterface;

  double radiusOfLaplacianForViscosity;
  double radiusOfLaplacianForViscosity_ratio;

  double radiusOfLaplacianForPressure;
  double radiusOfLaplacianForPressure_ratio;

  double radiusOfLaplacianForDiffusion;
  double radiusOfLaplacianForDiffusion_ratio;

  double maxRadius;
  double minRadius;

  double maxRadius_squared;
  double minRadius_squared;

  double maxRadius_ratio;
  double minRadius_ratio;

  double nZeroOfParticleNumberDensity; /* nZero(n0): standard particle number density */
  double nZeroOfGradient;
  double nZeroOfLaplacianForPressure;
  double nZeroOfLaplacianForViscosity;
  double nZeroOfLaplacianForDiffusion;
  int    nZeroOfRigidBody;
  double nZero_withConstantWeight;

  double nZeroOfnumberOfNeighborParticles;

  double lambdaOfLaplacianForPressure;
  double lambdaOfLaplacianForViscosity;
  double lambdaTimesNZeroForViscosity;
  double lambdaOfLaplacianForDiffusion;

  double relaxationCoefficientOfLambda;

  int    wallType;
  int    dummyWallType;
  int    typeNumberOfRigidParticle_forForcedMotion;

  double courantNumber;
  double diffusionNumber;


  double collisionDistance;
  double collisionDistance_ratio;
  double collisionDistance_squared;
  double collisionCoefficient;

  double thresholdOfParticleNumberDensity_ratio;

  int    numberOfParticleTypes;

  int    maxIterationNumberInIterationSolver;
  int    minIterationNumberInIterationSolver;

  double smallNumberForCheckingConvergenceInIterationSolver;

  double marginRatioForSettingCapacityOfLargeNeighborTable;
  double marginRatioForSettingCapacityOfSmallNeighborTable;

  int    capacityOfNeighborTable_large;
  int    capacityOfNeighborTable_small;
  int    maxArraySizeOfOccupiedNeighborTable_large;
  int    maxArraySizeOfOccupiedNeighborTable_small;

  double artificialViscosityCoefficient;


  double increaseRateOfDiagonalTermInCoefficientMatrixOfPressure;

  int    flagOfForcedMotionOfRigidBody;
  int    flagOfViscosityCalculation;

  int    flagOfHighViscosityCalculation;


  int    flagOfAutoSettingOfCapacityOfNeighborTable;
  int    flagOfExponentialExpressionInProfFile;

  int    flagOfCompressionOfProfFile;


  int    flagOfDivisionOfProfFile;
  int    flagOfBiCG;


  int    flagOfDivisionOfVtkFile;



  int    flagOfOutputOfTorqueFile;
  double numberOfRotations;


  int    flagOfOutputOfPressureFile;
  int    flagOfOutputOfAllWallParticle;
  int    flagOfDivisionOfPressureFile;

  int    numberOfPressureParticles;
  int    *numberListOfParticlesForOutputingPressure;

  int    upperLimitOfNumberOfProfFiles;
  int    upperLimitOfNumberOfPressureFiles;


  char   nameOfDataFile[SIZE_OF_FILE_NAME];
  char   nameOfGridFile[SIZE_OF_FILE_NAME];
  char   nameOfProfFile[SIZE_OF_FILE_NAME];

  char   nameOfVtkFile[SIZE_OF_FILE_NAME];


  char   nameOfDividedProfFile[SIZE_OF_FILE_NAME];

  char   nameOfDividedVtkFile[SIZE_OF_FILE_NAME];

  char   nameOfLogFile[SIZE_OF_FILE_NAME];
  char   nameOfSamplingDataFileForForcedMotionOfRigidBody[SIZE_OF_FILE_NAME];
  char   nameOfSamplingDataFileForForcedMotionOfWall[SIZE_OF_FILE_NAME];

  char   nameOfDesignationFileForWritingPressureFile[SIZE_OF_FILE_NAME];
  char   nameOfOutputPressureFile[SIZE_OF_FILE_NAME];
  char   nameOfOutputPressureFile_divided[SIZE_OF_FILE_NAME];


  char   nameOfDesignationFileForWritingTorqueFile[SIZE_OF_FILE_NAME];
  char   nameOfOutputTorqueFile[SIZE_OF_FILE_NAME];
  char   nameOfOutputTorqueFile_divided[SIZE_OF_FILE_NAME];

  char   nameOfMixingIndexFile[SIZE_OF_FILE_NAME];

  int    numberOfFluidParticles;
  int    numberOfRigidParticles;
  int    numberOfWallParticles;
  int    numberOfDummyWallParticles;
  int    numberOfGhostParticles;
  int    numberOfLiquidParticles[3];

  int    numberOfDesignatedParticles;

  int    numberOfNeighborTables;


  double relaxationCoefficientForSorSolver;

  double betaZeroOfParticles;

  int    flagOfBubbleCalculation;
  int    numperOfParticleForCalculatingBeta;
  char   nameOfBubbleInputFile[SIZE_OF_FILE_NAME];

  int    flagOfTanakaAndMasunagaModel;
  double valueOfGamma;
  double valueOfC;

  int    flagOfConcentrationCalculation;
  double DiffusionCoefficient;
  char   nameOfConcentrationInputFile[SIZE_OF_FILE_NAME];

  double averageConcentration[3];
  double maxConcentrationDifference[3];
  double mixingIndex[3];

  double pipeLength;

  double interfaceArea[10][10];  //for each particle type with each particle type
  double totalInterfaceArea[6];
  double totalLiquidInterfaceArea[6];

  int    flagOfInterfaceCalculation;
  int    flagOfInterfaceCorrection;
  int    flagOfInterfaceCompare;
  int    flagOfPeriodicX;

  int    flagOfSingleVortex;
  double alphaZero;

  int    threads;

  int    flagOfHideOutput;

  int    flagOfFabs; //Remove Again later, just for a test in Interface Calculation
  int    orderOfLsmps;
  int    flagOfCombinedInterfaceCalculation;

}structParameter;




typedef struct{

  int *list;
  int  count;

} structBucket;



typedef struct{

  double upperLimit[3];
  double upperLimit_previous[3];

  double lowerLimit[3];
  double lowerLimit_previous[3];

  double upperMarginRatio[3];
  double lowerMarginRatio[3];


  structBucket  ***bucket;

  int    capacityOfBucket;
  double marginRatioForSettingBucketCapacity;
  int    maxNumberOfParticlesInBucket;

  double bucketWidth[3];
  double bucketWidth_ratio[3];
  int    numberOfBuckets[3];

  int    flagOfAutoSettingOfDomainSize;
  int    flagOfAutoSettingOfBucketCapacity;
  int    flagOfOptimizationOfBucketMemorySize;


} structDomain;





typedef struct{

  int    numberOfComponents;
  int    *numberListOfComponents;
  double **relativePosition;
  double **relativePosition_rotated;
  double rotationalPosition[3];

  double rotationalMatrix[3][3];

  double translationalVelocity[3];
  double angularVelocity[3];

  double centerOfGravity[3];
  double centerOfGravity_previous[3];

  double centerOfWidth[3];

  double centerOfRotation[3];
  double centerOfRotation_previous[3];
  double centerOfRotation_initial[3];
  double centerOfRotation_directInput[3];

  int    particleType;
  int    howToSetCenterOfRotation;

  double quaternion[4];


  double maxPositionOfComponents[3];
  double minPositionOfComponents[3];

  int    *surfaceFlag;

} structObject;






typedef struct{

  int    numberOfSamplingData;

  double *samplingTime;

  double **translationalPosition;
  double **rotationalPosition;

  double interpolatedPosition_translational[3];
  double interpolatedPosition_rotational[3];

  double initialPosition_translational[3];
  double initialPosition_rotational[3];

  double startingTimeInSamplingData;
  double correctedSimulationTime;

  structObject object;

} structForcedMotion;



typedef struct{

  double real;
  double imaginary;

}structComplex;


structParticle          particle;
structDomain            domain;
structParameter         parameter;
structTimer             timer;
structPhysicalProperty  physicalProperty;
structObject            objectOfRigidBody;
structObject            objectOfWall;
structForcedMotion      forcedMotionOfRigidBody;
structForcedMotion      forcedMotionOfWall;


