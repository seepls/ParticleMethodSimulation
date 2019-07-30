
void
VISCOSITY_calculateViscosity( void );

void
VISCOSITY_setCoefficientMatrixOfViscosity(int iDim);

void
VISCOSITY_setSourceTerm( int iDim );

void
VISCOSITY_setTemporaryVelocityToZero( int iDim );

void
VISCOSITY_solveSimultaniousEquations( int iDim );

double
VISCOSITY_returnSmallNumberForCheckingConvergence( void );

void
VISCOSITY_setInitialSolution(void);
