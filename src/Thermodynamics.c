#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Thermodynamics.h"
#include "ConservedVariables.h"

#define GAMMA 1.4

void updatePressure(ConservedVariables* W)
{
    for (int i=0; i<W->nCells_i; i++)
    {
        for (int j=0; j<W->nCells_j; j++)
        {
            W->p[i][j] = (GAMMA-1)*(W->rhoE[i][j] - 1.0/(2.0*W->rho[i][j])*(pow(W->rhou[i][j],2)+pow(W->rhov[i][j], 2)));
        }
    }
}

