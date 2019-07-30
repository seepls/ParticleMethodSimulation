
void
TIMER_initializeTimeFunctions( void );

void
TIMER_setDtAutomatically( void );


void
TIMER_checkThatDtIsNotTooSmall( void );



double
TIMER_calculateDtWithCourantNumberInMind( void );


double
TIMER_calculateDtWithDiffusionNumberInMind( void );

void
TIMER_displayStateOfTimeStep_atAppropriateTime( void );


int
TIMER_checkWhetherItIsTimeToDisplayStateOfTimeStep( void );


void
TIMER_displayStateOfTimeStep( FILE *fp );


void
TIMER_displaySimulationTime( FILE *fp );


void
TIMER_displayCalculationTime( FILE *fp );

void
TIMER_displayRealTime( FILE *fp );


void
TIMER_calculateCalculationTime( void );


void
TIMER_putTimeForwardByDt( void );


int
TIMER_checkWhetherItIsTimeToWriteProfFile( void );

int
TIMER_checkWhetherItIsTimeToFinishProgram( void );

double
TIMER_getWTime(void);
