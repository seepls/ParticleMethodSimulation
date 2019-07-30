void
OTHER_endProgram( char *errorMessage );


void
OTHER_displayErrorMessage( char *errorMessage);


void
OTHER_normalEnd( void );


void
OTHER_changeParticleTypeIntoGhost( int iParticle );


void
OTHER_displayThatParticleBecameGhost( int iParticle );


void
OTHER_checkThatParticlesAreNotAllGhost( void );


void
OTHER_fillOneDimDoubleArrayWithZero( int num, double *array );


void
OTHER_fill2dimDoubleArrayWithZero( int nRow, int nColumn, double **array );


int
OTHER_isnan( double value );


int
OTHER_checkWhetherTheTwoPositionsAreNearyEqual( double *position1, double *position2 );
