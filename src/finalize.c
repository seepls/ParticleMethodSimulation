#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "file.h"
#include "copy.h"
#include "bucket.h"
#include "domain.h"
#include "neigh.h"
#include "density.h"
#include "timer.h"
#include "other.h"
#include "gravity.h"
#include "convection.h"
#include "collision.h"
#include "pressure.h"
#include "gradient.h"
#include "viscosity.h"
#include "sort.h"
#include "object.h"
#include "forcedMotion.h"
#include "memory.h"
#include "nzero.h"
#include "init.h"
#include "finalize.h"



void
FINALIZE_finalizeProgram( void ){


  NEIGH_displayRecommendedCapacityOfNeighborTable();
  BUCKET_displayRecommendedCapacityOfBucket();

  TIMER_displayCalculationTime( stderr );
  TIMER_displayCalculationTime( FpForLog );

  OTHER_normalEnd();

}


