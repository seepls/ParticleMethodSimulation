int*
MEMORY_allocateMemoryFor1dimIntArray( int num, char *variableName );


double*
MEMORY_allocateMemoryFor1dimDoubleArray( int num, char *variableName );
    

char*
MEMORY_allocateMemoryFor1dimCharArray( int num, char *variableName );


int**
MEMORY_allocateMemoryFor2dimIntArray( int nRow, int nColumn, char *variableName );


double**
MEMORY_allocateMemoryFor2dimDoubleArray( int nRow, int nColumn, char *variableName );


void
MEMORY_freeMemoryOf2dimIntArray( int **array, int nLine );


void
MEMORY_freeMemoryOf2dimDoubleArray( double **array, int nLine );


void
MEMORY_allocateMemoryForParticleStructure( void );


void
MEMORY_allocateMemoryForCoefficientMatrixOfPressure( void );
