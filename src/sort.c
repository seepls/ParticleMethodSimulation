#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "other.h"
#include "sort.h"



void
SORT_sortNeighborTable(  int  totalNumber
						 ,int  *numberOfNeighborParticles
						 ,int  **neighborTable
						 ,int  order ){

  int iParticle;

  if( order == ASCENDING_ORDER ){

    for(iParticle=0; iParticle < totalNumber; iParticle++){
      SORT_quickSort_inAscendingOrder(  iParticle, neighborTable, 0, numberOfNeighborParticles[iParticle]-1 );
    }

  }else if( order == DESCENDING_ORDER ){

    for(iParticle=0; iParticle < totalNumber; iParticle++){
      SORT_quickSort_inDescendingOrder(  iParticle, neighborTable, 0, numberOfNeighborParticles[iParticle]-1 );
    }

  }else{
    OTHER_endProgram("The last argument is not adequate. \nOrder must be ASCENDING_ORDER or DESCENDING_ORDER.\n[in function of SORT_sortNeighborTable()]");
  }

}



void
SORT_quickSort_inAscendingOrder(  int iParticle
								  ,int **neighborTable
								  ,int left
								  ,int right ){


  int  i,j;
  int  center;
  int  tmp;


  if(left < right){
    center = neighborTable[iParticle][(left+right)/2];

    i = left - 1; j = right + 1;


    while (1) {

      while (neighborTable[iParticle][++i] < center);
      while (neighborTable[iParticle][--j] > center);


      if ( i >= j ) break;

      tmp = neighborTable[iParticle][i]; 
      neighborTable[iParticle][i] = neighborTable[iParticle][j]; 
      neighborTable[iParticle][j] = tmp; 

    }

    SORT_quickSort_inAscendingOrder( iParticle, neighborTable, left, i-1);
    SORT_quickSort_inAscendingOrder( iParticle, neighborTable, j+1, right);

  }
}




void
SORT_quickSort_inDescendingOrder( int iParticle
								  ,int **neighborTable
								  ,int left
								  ,int right ){


  int i,j;
  int center;
  int tmp;


  if(left < right){
    center = neighborTable[iParticle][(left+right)/2];

    i = left - 1; j = right + 1;


    while (1) {

      while (neighborTable[iParticle][++i] > center);
      while (neighborTable[iParticle][--j] < center);


      if ( i >= j ) break;

      tmp = neighborTable[iParticle][i]; 
      neighborTable[iParticle][i] = neighborTable[iParticle][j]; 
      neighborTable[iParticle][j] = tmp; 

    }

    SORT_quickSort_inDescendingOrder( iParticle, neighborTable, left, i-1);
    SORT_quickSort_inDescendingOrder( iParticle, neighborTable, j+1, right);

  }
}
