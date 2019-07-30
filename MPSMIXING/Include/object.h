void
OBJECT_initializeObject( structObject *object, int particleType );


void
OBJECT_countNumberOfComponents( structObject *object );


int 
OBJECT_checkWhetherTheParticleIsComponent( int iParticle , structObject *object );


void
OBJECT_allocateMemoryForObjectStructure( structObject *object );


void
OBJECT_setNumberListOfComponents( structObject *object );


void
OBJCECT_addTheParticleInListOfComponents( int iParticle, structObject *object, int iComponent );


void
OBJECT_calculateCenterOfGravity( structObject *object );


void
OBJECT_calculateCenterOfWidth( structObject *object );


void
OBJECT_setInitialMaxMinPositionOfComponents(  double *maxPositionOfComponents
											 ,double *minPositionOfComponents
											  ,structObject   *object );


void
OBJECT_setCenterOfRotation( structObject *object );


void
OBJECT_memorizeInitialCenterOfRotation( structObject *object );


void
OBJECT_displayStateOfObject( structObject *object );


void
OBJECT_setRelativePosition( structObject *object );


void
OBJECT_reshapeObject( structObject *object );


void
OBJECT_initializeRotatedRelativePosition(
										 double       *rotatedRelativePosition
										 ,int          iComponent
										 ,structObject *object
										 );


void
OBJECT_rotateRelativePosition(
							  double       *rotatedRelativePosition
							  ,structObject *object
							  );


void
OBJECT_setVelocityOfObject( structObject *object );


void
OBJECT_updateObjectProperty( structObject *object );


void
OBJECT_resetMatrix( double matrix[3][3] );


void
OBJECT_writeObjectPositionInFile( structObject *object );


void
OBJECT_calculateRotationalPosition( double *rotationalPosition, double *quaternion);
