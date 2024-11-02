#ifndef CONSERVED_VARIABLES_H
#define CONSERVED_VARIABLES_H

#include "mesh.h"

typedef struct
{
    int n_ghost_layers;
    int nCells_i; // number of cells in the i direction including ghost cells
    int nCells_j; // number of cells in the j direction including ghost cells
    float** rho;
    float** rhou;
    float** rhov;
    float ** rhoE;
    float ** p;
}ConservedVariables;

ConservedVariables* createConservedVariables(MESH* mesh, int n_ghost_layers,float alpha, float rho_inf, float p_inf, float M_inf);
void FreeConservedVariables(ConservedVariables* W);
void initConservedVariables(ConservedVariables* W, float alpha, float p_inf, float rho_inf, float M_inf);


#endif