void
INTERFACE_CalculateInterfaceArea(void);

void
INTERFACE_CalculateInterfaceAreaBetweenIAndJB(int iType, int jType);

void
INTERFACE_CalculateInterfaceAreaBetweenIAndJ(int iType, int jType);

int
INTERFACE_setFabs(int iParticle, int numberOfNeighborParticles, int *neighborTable, int *neighborTablePeriodic, int iType, int jType);

void
INTERFACE_calculateInterfaceCompare(void);

void 
INTERFACE_CalculateInterfaceAreaBetweenIAndJLSMPS(int iType, int jType);

void 
INTERFACE_CalculateInterfaceAreaBetweenIAndJSplit(int iType, int jType);

void
//INTERFACE_calculateCorrectionMatrix(int iParticle, int numberOfNeighborParticles, int *neighborTable, int *neighborTablePeriodic, double **corrMatrixIJ);
INTERFACE_calculateCorrectionMatrix(int iParticle, int numberOfNeighborParticles, int *neighborTable, int *neighborTablePeriodic, double corrMatrixIJ[MAX_DIM][MAX_DIM]);

void
INTERFACE_matrixInverse2x2b(double inverse[10][10]);

void
INTERFACE_matrixInverse2x2(double matrix[MAX_DIM][MAX_DIM], double inverse[MAX_DIM][MAX_DIM]);

void
INTERFACE_matrixInverse3x3(double matrix[MAX_DIM][MAX_DIM], double inverse[MAX_DIM][MAX_DIM]);

void
INTERFACE_initializeInterfaceFile(void);

void
INTERFACE_initializeInterfaceFileCompare(void);

void
INTERFACE_writeInterfaceFile(int iType, int jType);

void
INTERFACE_writeInterfaceFileCompare(int iType, int jType);

void
INTERFACE_Threshold(double nablaF[MAX_DIM], double lowerThreshold, double upperThreshold, double criteria);
