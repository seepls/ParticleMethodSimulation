void
NEIGH_initializeNeighborTable( void );


void
NEIGH_setCapacityOfNeighborTableAutomatically( void );


void
NEIGH_setNeighborTable( double **position );

void
NEIGH_setNeighborTablePeriodic( double **position );


void
NEIGH_recordMaxArraySizeOfOccupiedNeighborTable( void );



void
NEIGH_checkCapacityOfNeiborTable( 
								 int  iParticle
								 ,int  *numberOfNeighborParticles
								 ,int  **neighborTable
								 ,int  capacity
								 ,char *nameOfTable
								 );

void
NEIGH_addParticleInNeighborTable( 
								 int iParticle
								 ,int jParticle
								 ,int *numberOfNeighborParticles
								 ,int **neighborTable
								 );


void
NEIGH_resetNeighborTable( void );


void
NEIGH_selectNeighborTable(  int **numberOfNeighborParticles
							,int ***neighborTable
							,double radius_ratio
							,int ***neighborTablePeriodic
							);



void
NEIGH_displayCapacityOfNeighborTable( void );


void
NEIGH_displayErrorMessageForMemoryOverFlow(
										   int  iParticle
										   ,int  *numberOfNeighborParticles
										   ,int  **neighborTable
										   ,char *nameOfTable
										   );


void
NEIGH_displayRecommendedCapacityOfNeighborTable( void );


void
NEIGH_allocateMemoryForNeighborTable( void );
