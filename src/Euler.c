#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "mesh.h"
#include "writer.h"
#include "ConservedVariables.h"
#include "Thermodynamics.h"
#include "BoundaryCondition.h"

int main() {

    // Parameters
    char* meshPath = "NACA0012grids/65x65.x";
    char* outputFile = "output/Flow.dat";
    int n_ghost_layers = 2;

    float rho_inf = 1.0;
    float p_inf = 1.0;
    float M_inf = 0.3;
    float alpha = 5*PI/180.0;

    // 1. Read the mesh and compute geometrical quantities (normal, face lengths, cell volumes)
    MESH* MESH = LoadMESH(meshPath);

    // 2. Create and initialize the conserved variables
    ConservedVariables* W = createConservedVariables(MESH, n_ghost_layers, alpha, rho_inf, p_inf, M_inf);


    // 3. Create boundary conditions

    // Remark :  index i and j are given in the inner mesh framework -> not considering the ghost cells
    // -1 is used to say that there is no condition in the given direction
    BoundaryCondition* bc_wall = createBoundaryCondition(WALL, I_DIRECTION, 0, -1);
    BoundaryCondition* bc_farfield = createBoundaryCondition(FARFIELD, I_DIRECTION, MESH->ni-2, -1);
    BoundaryCondition* bc_connect1 = createBoundaryCondition(CONNECT, J_DIRECTION, -1, 0);
    BoundaryCondition* bc_connect2 = createBoundaryCondition(CONNECT, J_DIRECTION, -1, MESH->nj-2);
    // -2 for the max index since we work on cells not on vertices

    applyWallBoundary(W, MESH, bc_wall); // test of bcwall
    updatePressure(W);


    // 4. Write results
    WriteTecplotFile(outputFile, MESH, W);

    // Free the allocated arrays at the end 
    FreeMESH(MESH);
    FreeConservedVariables(W);
    free(bc_connect1);
    free(bc_connect2);
    free(bc_wall);
    free(bc_farfield);
    return 0; // Return success
}