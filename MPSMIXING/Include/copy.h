
void
COPY_setInitialCoodinatesAndVelocities( void );

void
COPY_storeInitialParticleProperty( void );


void
COPY_updateParticleProperty( void );


void
COPY_copy2dimDoubleArray( int, int, double **, double ** );

void
COPY_copy1dimDoubleArray( int, double *, double * );


void
COPY_copy1dimIntArray( int num, int *overwrittenArray, int *storedArray );

void
COPY_copy1dimDoubleArrayfrom2dimDoubleArray( int nRow, int nColumn, double *overwrittenArray, double **storedArray);

void
COPY_copy2dimDoubleArrayfrom1dimDoubleArray( int nRow, int nColumn, double **overwrittenArray, double *storedArray);
