
FILE*
FILE_openFile(char *filename, char *mode);


void
FILE_closeFile( FILE *fp, char *fileName);


void
FILE_readInputFiles( int argumentCount,char **argumentVector );


void
FILE_openLogFile( void );


void
FILE_setNameOfInputFile( int argumentCount,char **argumentVector );


void
FILE_setNameOfOutputFile( int argumentCount,char **argumentVector );


void
FILE_displayArgumentOfCommandLine( int argumentCount, char **argumentVector );


void
FILE_setDefaultNameOfInputFile( void );


void
FILE_setDefaultNameOfOutputFile( void );


void
FILE_readGridFile( void );


void
FILE_readDesignationFileForWritingPressureFile( void );


void
FILE_countNumberOfParticlesEachType( void );


void
FILE_displayNumberOfParticles( void );


void
FILE_scanDouble( FILE *fp, double *scanedValue, char *parameterName);


void
FILE_scanInt( FILE *fp, int *scanedValue, char *parameterName);


void
FILE_scanChar( FILE *fp, char *scanedValue, char *parameterName);


void
FILE_scanOnOff(FILE *fp, int *flag, char *variable_name);


void
FILE_scanTypeOfOutputInterval(FILE *fp, int *flag, char *variable_name);


void
FILE_scanMotionTypeOfRigidBody(FILE *fp, int *flag, char *variable_name);


void
FILE_scanDirectionOfWavePropagation(FILE *fp, int *flag, char *variable_name);


void
FILE_scanKindOfSolverForSimultaneousEquation( FILE *fp);


char*
FILE_returnOnOff(int flag);


char*
FILE_returnKindOfSolverForSimultaniousEquation( void );


char*
FILE_returnDim(int iDim);



char*
FILE_returnTypeOfOutputInterval(int flag);



char*
FILE_returnMotionTypeOfRigidBody(int flag);


char*
FILE_returnDirectionOfWavePropagation(int flag);


void
FILE_skipComment( FILE *fp );


void
FILE_checkEndOfFile(int status, char *name);


void
FILE_displayGridInformation( void );


void
FILE_writeCalculationResultInFile( void );


void
FILE_writeProfFile( void );


void
FILE_writePressureFile( void );


void
FILE_outputPressureOfAllWallParticles( void );


void
FILE_outputPressureOfOnlyDesignatedParticles( void );

void
FILE_writeTorqueFile( void );


void
FILE_outputTorqueOfAllRigidParticles( void );


void
FILE_outputTorqueOfOnlyDesignatedParticles( void );


FILE*
FILE_openProfFile( void );


FILE*
FILE_openPressureFile( void );


void
FILE_compressProfFile( void );


void
FILE_readDataFile( void );


void
FILE_displayReadDataFile( void );


FILE*
FILE_openVtkFile( void );

void
FILE_makeVtkFile( void );
