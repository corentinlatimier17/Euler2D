#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "BoundaryCondition.h"

BoundaryCondition* createBoundaryCondition(BoundaryConditionType type, BoundaryDirection direction, int index_i, int index_j)
{
    BoundaryCondition* bc = (BoundaryCondition*)malloc(sizeof(BoundaryCondition));
    bc->direction = direction;
    bc->type = type;
    bc->index_i = index_i;
    bc->index_j = index_j;
}

void applyWallBoundary(ConservedVariables* W, MESH* mesh, BoundaryCondition* bc_wall) 
{
    // no slip wall bc
    if (bc_wall->direction == I_DIRECTION) // wall is at constant i
    {
        int i_wallW = bc_wall->index_i + W->n_ghost_layers; // i_wallW is the index in W grids
        int i_wall = bc_wall->index_i; // i_wall is the index in inner mesh (without ghost cells)

        for (int j=0; j<mesh->nj-1; j++) // Loop over j direction (loop over cells)
        {
            int jw = j + W->n_ghost_layers;
            Normal2D normal = getNormal(i_wall, j, i_wall, j+1, mesh);
            if (i_wall==mesh->ni-2) // we get -n, not n
            {
                normal.nx= -1.0*normal.nx;
                normal.ny = -1.0*normal.ny;
            }
            float V2n = W->rhou[i_wallW][jw]/W->rho[i_wallW][jw]*normal.nx + W->rhov[i_wallW][jw]/W->rho[i_wallW][jw]*normal.ny; 
            
            float ubc = W->rhou[i_wallW][jw]/W->rho[i_wallW][jw] -2.0*V2n*normal.nx;
            float vbc = W->rhov[i_wallW][jw]/W->rho[i_wallW][jw] -2.0*V2n*normal.ny;

            if (i_wall == 0)
            {
                W->rhou[i_wallW-1][jw]=W->rho[i_wallW][jw]*(ubc);
                W->rhov[i_wallW-1][jw]=W->rho[i_wallW][jw]*(vbc);
                W->rho[i_wallW-1][jw] = W->rho[i_wallW][jw];
                W->rhoE[i_wallW-1][jw] = W->rhoE[i_wallW][jw];

                W->rhou[i_wallW-2][jw]=W->rhou[i_wallW-1][jw];
                W->rhov[i_wallW-2][jw]=W->rhov[i_wallW-1][jw];
                W->rhoE[i_wallW-2][jw]=W->rhoE[i_wallW-1][jw];
                W->rho[i_wallW-2][jw]=W->rho[i_wallW-1][jw];
            }
            else if (i_wall == mesh->ni-2) // -2 since we consider a cell
            {
                    W->rhou[i_wallW + 1][jw] = W->rho[i_wallW][jw] * ubc;
                    W->rhov[i_wallW + 1][jw] = W->rho[i_wallW][jw] * vbc;
                    W->rho[i_wallW + 1][jw] = W->rho[i_wallW][jw];
                    W->rhoE[i_wallW + 1][jw] = W->rhoE[i_wallW][jw];

                    W->rhou[i_wallW + 2][jw] = W->rhou[i_wallW + 1][jw];
                    W->rhov[i_wallW + 2][jw] = W->rhov[i_wallW + 1][jw];
                    W->rhoE[i_wallW + 2][jw] = W->rhoE[i_wallW + 1][jw];
                    W->rho[i_wallW + 2][jw] = W->rho[i_wallW + 1][jw];
            } 
        }
    }
    if (bc_wall->direction == J_DIRECTION) // wall is at constant j
    {
        int j_wallW = bc_wall->index_j + W->n_ghost_layers; // i_wallW is the index in W grids
        int j_wall = bc_wall->index_j; // i_wall is the index in inner mesh (without ghost cells)

        for (int i=0; i<mesh->ni-1; i++) // Loop over i direction
        {
            int iw = i + W->n_ghost_layers; // index in W grids
            Normal2D normal = getNormal(i, j_wall, i+1, j_wall, mesh);
            if (j_wall==0) // we get -n, not n
            {
                normal.nx= -1.0*normal.nx;
                normal.ny = -1.0*normal.ny;
            }

            float V2n = W->rhou[iw][j_wallW]/W->rho[iw][j_wallW]*normal.nx + W->rhov[iw][j_wallW]/W->rho[iw][j_wallW]*normal.ny; 
            
            float ubc = W->rhou[iw][j_wallW]/W->rho[iw][j_wallW] -2.0*V2n*normal.nx;
            float vbc = W->rhov[iw][j_wallW]/W->rho[iw][j_wallW] -2.0*V2n*normal.ny;

            if (j_wall == 0)
            {
                W->rhou[iw][j_wallW-1]=W->rhou[iw][j_wallW]*(ubc);
                W->rhov[iw][j_wallW-1]=W->rhov[iw][j_wallW]*(vbc);
                W->rho[iw][j_wallW-1] = W->rho[iw][j_wallW];
                W->rhoE[iw][j_wallW-1] = W->rhoE[iw][j_wallW];

                W->rhou[iw][j_wallW-2]=W->rhou[iw][j_wallW-1];
                W->rhov[iw][j_wallW-2]=W->rhov[iw][j_wallW-1];
                W->rhoE[iw][j_wallW-2]=W->rhoE[iw][j_wallW-1];
                W->rho[iw][j_wallW-2]=W->rho[iw][j_wallW-1];
            }
            else if (j_wall == mesh->nj-2) // -2 since we consider a cell
            {
                W->rhou[iw][j_wallW+1]=W->rhou[iw][j_wallW]*(ubc);
                W->rhov[iw][j_wallW+1]=W->rhov[iw][j_wallW]*(vbc);
                W->rho[iw][j_wallW+1] = W->rho[iw][j_wallW];
                W->rhoE[iw][j_wallW+1] = W->rhoE[iw][j_wallW];

                W->rhou[iw][j_wallW+2]=W->rhou[iw][j_wallW+1];
                W->rhov[iw][j_wallW+2]=W->rhov[iw][j_wallW+1];
                W->rhoE[iw][j_wallW+2]=W->rhoE[iw][j_wallW+1];
                W->rho[iw][j_wallW+2]=W->rho[iw][j_wallW+1];
            }
        }
    }
}