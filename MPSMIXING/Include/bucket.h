void
BUCKET_initializeBucket( void );


void
BUCKET_storeParticlesInBuckets( double **position );


void
BUCKET_countTheNumberOfParticleInEachBucket( double **position );


void
BUCKET_allocateMemoryForParticleListOfBucket( double **position );


void
BUCKET_freeMemoryOfParticleList( void );


void
BUCKET_storeParticlesInBucketsNORMAL( double **position );


void
BUCKET_setWidthOfBucket( void );


void
BUCKET_setNumberOfBuckets( void );


void
BUCKET_setCapacityOfBucketAutomatically( void );


void
BUCKET_displayInformationOfBucket( void );


void
BUCKET_allocateMemoryForBuckets( void );


void
BUCKET_checkCapacityOfTheBucket( int currentNumberOfStoredParticles, int iX, int iY, int iZ );


void
BUCKET_displayErrorMessageForMemoryOverFlow( int iX, int iY, int iZ );


void
BUCKET_storeTheParticleInBucket(int iParticle, int endOfList, int iX, int iY, int iZ);


void
BUCKET_resetParticleCounterToZero( void );


void 
BUCKET_findBucketWhereParticleIsStored( 
									   int *iX, int *iY, int *iZ
									   ,int            iParticle
									   ,double        **position
									   );


void
BUCKET_checkErrorOfStorePlace(int iX, int iY, int iZ, int iParticle, double **position);


void
BUCKET_freeBuckets( void );


void
BUCKET_recordMaxNumberOfParticlesInBucket( void );


void
BUCKET_displayRecommendedCapacityOfBucket( void );
