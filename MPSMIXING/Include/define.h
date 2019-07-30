#define SIZE_OF_FILE_NAME 1024
#define MAX_DIM        3

#define TRUE           1
#define FALSE          0

#define EXCEPTION          -1
#define EXCEPTION_OCCURED  -1
#define YES            TRUE
#define NO             FALSE

#define ON             TRUE
#define OFF            FALSE

#define XDIM           0
#define YDIM           1
#define ZDIM           2

#define GHOST          -1
#define GHOST_OR_DUMMY -1

#define SURFACE_PARTICLE 1
#define INNER_PARTICLE   0

#define CONNECTED      1
#define NOT_CONNECTED  2
#define CHECKED        3

#define IN_DOMAIN      1
#define OUT_OF_DOMAIN  2

#define ASCENDING_ORDER  1
#define DESCENDING_ORDER 2

#define TOO_SMALL_VALUE 1.0E-9

#define     CG_SOLVER 1
#define   BICG_SOLVER 2
#define    SOR_SOLVER 3
#define JACOBI_SOLVER 4



#define NAME_OF_DATA_FILE         "input.data"
#define NAME_OF_GRID_FILE         "input.grid"
#define NAME_OF_PROF_FILE         "output.prof"
#define NAME_OF_DIVIDED_PROF_FILE "output_"
#define NAME_OF_LOG_FILE          "output.log"

#define NAME_OF_BUBBLE_FILE     "input.bubble"
#define NAME_OF_CONCENTRATION_FILE  "input.concentration"

#define NAME_OF_VTK_FILE          "output.vtk"
#define NAME_OF_DIVIDED_VTK_FILE  "output_"

#define NAME_OF_TORQUE_FILE          "output.torque"
#define NAME_OF_DIVIDED_TORQUE_FILE  "output_"

#define NAME_OF_MIXING_INDEX_FILE    "mixing.index"



#define SIMULATION_TIME 1
#define TIME_STEP       2

#define THRESHOLD_CALCULATION_TIME  1.0

#ifndef M_PI
#define	M_PI	3.1415926535897932384626433832795
#endif


#define FREE_MOTION       1
#define FORCED_MOTION     2

#define FROM_MINUS_TO_PLUS 1
#define FROM_PLUS_TO_MINUS 2

#define CENTER_OF_GRAVITY  1
#define CENTER_OF_WIDTH    2
#define DIRECT_INPUT       3


#define FLUID_REGION       1
#define WALL_REGION        2
#define INSIDE_OF_TANK     3
#define OUTSIDE_OF_TANK    4
#define GHOST_REGION       -1


#define UPPER_THE_WAVE_SURFACE 1
#define UNDER_THE_WAVE_SURFACE 2

#define BETWEEN_PARTICLES     1
#define NOT_BETWEEN_PARTICLES 2


/*--- START___difinitions for freeMotion of Rigid body ---*/

#define SMALL_VALUE_FOR_CHECKING_CONVERGENCE_OF_SOLVER_FOR_EIGEN_VALUE   1.0E-20
#define MAX_ITERATION_NUMBER_IN_JACOBI_SOLVER_FOR_EIGEN_VALUE            1000

#define CONVERGENT                                                       1
#define NOT_CONVERGENT                                                   2

#define NORMAL_CONVERGENCE                                               1
#define ITERATION_NUMBER_EXCEDDED_THE_LIMIT                              2


#define FORCE_BY_AREA_INTEGRATION_OF_PRESSURE                            1
#define FORCE_BY_MOMENTUM_CHANGE                                         2

/*--- END___difinitions for freeMotion of Rigid body --- */

#define MIRROR_PARTICLE      TRUE
#define NOT_MIRROR_PARTICLE  FALSE
#define MIRROR_OF_NONE          0
#define MIRROR_OF_UPPER_X_PLANE 1
#define MIRROR_OF_LOWER_X_PLANE 2
#define MIRROR_OF_UPPER_Y_PLANE 3
#define MIRROR_OF_LOWER_Y_PLANE 4
#define MIRROR_OF_UPPER_Z_PLANE 5
#define MIRROR_OF_LOWER_Z_PLANE 6

 

