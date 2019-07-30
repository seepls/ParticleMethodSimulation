void
DOMAIN_initializeDomain( void );


void
DOMAIN_setRangeOfDomain( void );


int
DOMAIN_checkWhetherParticleIsInDomain( int iParticle, double **position );

void
DOMAIN_relocateParticlePeriodic( int iParticle, double **position );


void
DOMAIN_displayUpdatedDomainSize( void );
