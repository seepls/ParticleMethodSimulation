void
SOLVER_solveSimultaniousEquations( void );


double
SOLVER_returnSmallNumberForCheckingConvergence( void );


void
SOLVER_setInitialSolution( void );


int
SOLVER_conjugateGradientMethod( 
								int      totalNumber
								,double **matrix
								,double  *solution
								,double  *sourceTerm
								,int     *flagOfBoundaryCondition
								,int     *numberOfNeighborParticles
								,int     **neighborTable
								,int      maxIteration
								,int      minIteration
								,double   smallNumberForCheckingConvergence
								);


void
SOLVER_displayStateOfConvergence( FILE *fp, int iIteration, int minIteration );


void
SOLVER_displayWarnigMessageForExcessOfIterationNumber( int iIteration, int maxIteration );


double
SOLVER_calculateNormOfSorceTerm(  
								 int    totalNumber
								,double *sourceTerm
								,int    *flagOfBoundaryCondition
								 );


int
SOLVER_checkConvergence(  
						 int    totalNumber
						,int    *flagOfBoundaryCondition
						,double *remainder
						,double normOfSourceTerm
						,int    iIteration
						,int    minIteration
						,double smallNumberForCheckingConvergence
						 );

int
SOLVER_checkConvergence_NORMAL(  
							   int    totalNumber
							   ,int    *flagOfBoundaryCondition
							   ,double *remainder
							   ,int    iIteration
							   ,int    minIteration
							   ,double smallNumberForCheckingConvergence
							   );


int
SOLVER_findIndexOfIParticleInNeighborListOfJParticle( int iParticle,  int jParticle);

