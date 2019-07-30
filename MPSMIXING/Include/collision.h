void
COLLISION_calculateCollisionBetweenParticles(void);


void
COLLISION_displayThatCollisionOccured(
									  int
									  ,int
									  ,double
									  ,double
									  );



void
COLLISION_countCollision(
						 int *coutOfCollision
						 ,int *coutOfCollisionBetween_fluidAndFluid
						 ,int *coutOfCollisionBetween_fluidAndWall
						 ,int *coutOfCollisionBetween_fluidAndDummyWall
						 ,int  iParticle
						 ,int  jParticle
						 );


void
COLLISION_displayCountOfCollision( 
								  int coutOfCollision
								  ,int coutOfCollisionBetween_fluidAndFluid
								  ,int coutOfCollisionBetween_fluidAndWall
								  ,int coutOfCollisionBetween_fluidAndDummyWall
								  );


