#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ConservedVariables.h"
#include "mesh.h"
#include "Thermodynamics.h"


ConservedVariables* createConservedVariables(MESH* mesh, int n_ghost_layers, float alpha, float rho_inf, float p_inf, float M_inf)
{
    ConservedVariables* W = (ConservedVariables*)malloc(sizeof(ConservedVariables));
    // init n_ghost_layers attribute
    W->n_ghost_layers = n_ghost_layers;
    W->nCells_i = 2*n_ghost_layers + mesh->ni-1;
    W->nCells_j = 2*n_ghost_layers + mesh->nj-1;

    // Initialize the conserved variables grid with error handling
    W->rho = (float**)malloc(W->nCells_i * sizeof(float*));
    W->rhou = (float**)malloc(W->nCells_i * sizeof(float*));
    W->rhov = (float**)malloc(W->nCells_i * sizeof(float*));
    W->rhoE = (float**)malloc(W->nCells_i * sizeof(float*));
    W->p = (float**)malloc(W->nCells_i * sizeof(float*));

    if (W->rho == NULL || W->rhou == NULL || W->rhov == NULL || W->rhoE == NULL || W->p == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for one of the main arrays for the conserved variables\n");
        // Free any allocated memory for consistency
        free(W->rho);
        free(W->rhou);
        free(W->rhov);
        free(W->rhoE);
        free(W->p);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each grid and initialize to 0.0f
    for (int i = 0; i < W->nCells_i; i++) {
        W->rho[i] = (float*)malloc(W->nCells_j * sizeof(float));
        W->rhou[i] = (float*)malloc(W->nCells_j * sizeof(float));
        W->rhov[i] = (float*)malloc(W->nCells_j * sizeof(float));
        W->rhoE[i] = (float*)malloc(W->nCells_j * sizeof(float));
        W->p[i] = (float*)malloc(W->nCells_j * sizeof(float));

        // Check if any allocation failed for the current row
        if (W->rho[i] == NULL || W->rhou[i] == NULL || W->rhov[i] == NULL || W->rhoE[i] == NULL || W->p[i] == NULL) {
            fprintf(stderr, "Error: Memory allocation failed for one of the arrays of the conserved variables at row %d\n", i);

            // Free previously allocated rows for each array
            for (int k = 0; k < i; k++) {
                free(W->rho[k]);
                free(W->rhou[k]);
                free(W->rhov[k]);
                free(W->rhoE[k]);
                free(W->p[k]);
            }
            // Free the main arrays
            free(W->rho);
            free(W->rhou);
            free(W->rhov);
            free(W->rhoE);
            free(W->p);
            
            exit(EXIT_FAILURE);
        }
    }
    initConservedVariables(W, alpha, p_inf, rho_inf, M_inf);
    return W;
}

void initConservedVariables(ConservedVariables* W, float alpha, float p_inf, float rho_inf, float M_inf)
{
    for (int i=0; i<W->nCells_i; i++)
    {
        for (int j=0; j<W->nCells_j; j++)
        {
            W->rho[i][j] = rho_inf;
            W->p[i][j] = p_inf;

            float v_inf = sqrt(GAMMA*p_inf/rho_inf);
            W->rhou[i][j] = W->rho[i][j]*v_inf*cos(alpha);
            W->rhov[i][j] = W->rho[i][j]*v_inf*sin(alpha);

            W->rhoE[i][j] = W->p[i][j]/(GAMMA-1) + 1.0/(2.0*W->rho[i][j])*(pow(W->rhou[i][j],2)  + pow(W->rhov[i][j],2));
        }
    }
}



void FreeConservedVariables(ConservedVariables* W) {
    if (W == NULL) return;

    // Free each row of the 2D arrays in W
    for (int i = 0; i < W->nCells_i; i++) {
        if (W->rho[i] != NULL) {
            free(W->rho[i]);
        }
        if (W->rhou[i] != NULL) {
            free(W->rhou[i]);
        }
        if (W->rhov[i] != NULL) {
            free(W->rhov[i]);
        }
        if (W->rhoE[i] != NULL) {
            free(W->rhoE[i]);
        }
        if (W->p[i] != NULL) {
            free(W->p[i]);
        }

    }

    // Free the main arrays themselves
    free(W->rho);
    free(W->rhou);
    free(W->rhov);
    free(W->rhoE);
    free(W->p);

    // Free the ConservedVariables structure itself
    free(W);
}

