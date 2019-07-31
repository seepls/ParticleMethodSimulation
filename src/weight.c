#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "other.h"
#include "weight.h"


double
WEIGHT_calculateWeightFunction( double distance, double radius ){

  double weightIJ;

  if ( distance < radius){

      weightIJ = radius/distance - 1.0;

  }else{

    weightIJ = 0.0;

  }

  return(weightIJ);

}

int
WEIGHT_calculateWeightZeroorOneFunction( double distance, double radius ){
    
    int weightIJ;
    
    if ( distance < radius){
        
        weightIJ = 1;
        
    }else{
        
        weightIJ = 0;
        
    }
    
    return(weightIJ);
    
}


