void
NZERO_calculateNZeroAndLambda( double averageDistanceBetweenParticles );


void
NZERO_setCoordinateOfVirtualParticles(  
									  double **positionOfVirtualParticle
									  ,int    *particleNumberOfCenterParticle
									  ,int    *numberOfVirtualParticles
									  ,int    nX
									  ,int    nY
									  ,int    nZ
									  ,double averageDistanceBetweenParticles
									  );


void
NZERO_setZeroToNZeroAndLambda( void );


void
NZERO_setSizeOfVirtualDomain(  
							 int *nX
							 ,int *nY
							 ,int *nZ
							 );


double**
NZERO_allocateMemoryForVirtualParticles( int nX, int nY, int nZ );


void
NZERO_displayCalculatedNZeroAndLambda( void );
