void
PRESSURE_calculatePressure( void );


void
PRESSURE_setSourceTerm( void );


void
PRESSURE_setCoefficientMatrixOfPressure( void );


void
PRESSURE_setMinusPressureZero( void );


void
PRESSURE_setMinimumPressureAmongNeighbors( void );


int
PRESSURE_checkNecessityOfPressureCalculation( void );


void
PRESSURE_setDirichletBoundaryCondition( void );


void
PRESSURE_setDirichletBoundaryCondition_twice(void);


void
PRESSURE_findParticlesBeingOnFreeSurface( void );

void
PRESSURE_findParticlesBeingOnFreeSurfaceTheNumberOfNeighParticleBased( void );


void
PRESSURE_checkThatDirichletBoundaryConditionIsConnected( void );


void
PRESSURE_ifBoundaryConditionIsNotConnected_increaseDiagonalElementsOfCoefficientMatrixOfPressureToSolveSimultaneousEquationsWithoutBoundaryCondition( int *checkArray );


void
PRESSURE_displayWarnigMessageForNoDirichletBoundaryCondition( int iParticle );
/*******************************************10.22追加******************************************/
void
PRESSURE_correctPressure( void );
/*********************************************************************************************/
