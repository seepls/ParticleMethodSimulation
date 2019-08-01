#include <stdio.h>
#include <math.h>
#include <string.h>

#include "define.h"
#include "struct.h"
#include "weight.h"
#include "extern.h"
#include "bubble.h"
#include "distance.h"
#include "buoyancy.h"
#include "file.h"
#include "gravity.h"
#include "copy.h"

void
Buoyancy_calculateBuyancy( void ){
    
    int iParticle;
    int iDim;
    double voidrate;
    
    for(iParticle= 0; iParticle < particle.totalNumber; iParticle++){
        
        if ( particle.type[iParticle]==GHOST ) continue;
        if ( particle.type[iParticle]==parameter.dummyWallType ) continue;
        if ( particle.type[iParticle]==parameter.wallType ) continue;
        
        voidrate=(particle.moleOfBubbles[iParticle]*physicalProperty.gasConstant*physicalProperty.temperature)/(((particle.pressure[iParticle])+physicalProperty.headPressure)*(Bubble_calculateVolume(particle.averageDistance)));
        
        
        for(iDim= 0; iDim < NumberOfDimensions; iDim++){
            
            particle.velocity[iDim][iParticle] -= ((physicalProperty.massDensity[0]-physicalProperty.massDensityOfBubble)/physicalProperty.massDensity[0]) * voidrate * physicalProperty.gravity[iDim] * timer.dt;
        
        }
        
    }
    
}
