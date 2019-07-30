double
MATHE_average( double a1, double a2 );


void
MATHE_multiplyMatrixByVector(  
							 int    totalNum
							 ,double *answer
							 ,double **matrix
							 ,double *vector
							 ,int    *numberOfNeighborParticles
							 ,int    **neighborList
							 ,int    *flagOfBoundaryCondition
							 );


void
MATHE_multiplyMatrixByMatrix_matrixSizeIsThree( double answer[3][3]
												,double matrix_left[3][3]
												,double matrix_right[3][3]
												);



void
MATHE_addVectorToVector( int num, double *answer, double *vector1, double *vector2, double a);


void
MATHE_add2dimArrayTo2dimArray( int nRow, int nColumn, double **answer, double **array1, double **array2 );


void
MATHE_subtractVectorFromVector( int num, double *answer, double *vector1, double *vector2, double a);


double
MATHE_innerProductOfVectors( double *vector1, double *vector2, int num);


double 
MATHE_innerProductOfVectorsForIterationSolver( double *vector1, double *vector2, int num,  int *flagOfBoundaryCondition);


double
MATHE_linearInterpolation(  double value2, double value1,  double m2, double m1, double n2, double n1 );


double
MATHE_absoluteValueOfVector( double *vector, int nDimension );


void
MATHE_outerProduct( 
				   double *answer
				   ,double *vector1
				   ,double *vector2
				   );


void
MATHE_normalizeVector( double *vector, int nDimension );


void
MATHE_linearTransform( double *answer, double matrix[3][3], double *vector );


void
MATHE_matrixInverse2x2(double **matrix, double **inverse);
