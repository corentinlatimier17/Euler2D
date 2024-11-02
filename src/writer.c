#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"
#include "ConservedVariables.h"

void WriteTecplotFile(const char* filename, MESH* mesh, ConservedVariables* W) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    // Write Tecplot header
    fprintf(file, "Title = \"Flow solution\"\n");
    fprintf(file, "Variables = X, Y, CellVolume, I, J, rhou, rhov, rho, rhoE, p\n");
    fprintf(file, "ZONE T=\"BLOCK1\", I=%d, J=%d\n, DATAPACKING=POINT\n", mesh->ni, mesh->nj);

    int n_ghosts_layers = W->n_ghost_layers;

    // Write the mesh coordinates and variables
    for (int j = 0; j < mesh->nj; j++) { // loop over nodes, not over cells -> take into account the ghost cells for the W variables
        for (int i = 0; i < mesh->ni; i++) 
        {
            float CellVolume = 0.0;
            int iW = i+n_ghosts_layers;
            int jW = j+n_ghosts_layers;

            float rho = 0.25*(W->rho[iW][jW]+W->rho[iW][jW-1]+W->rho[iW-1][jW-1]+W->rho[iW-1][jW]);
            float rhou = 0.25*(W->rhou[iW][jW]+W->rhou[iW][jW-1]+W->rhou[iW-1][jW-1]+W->rhou[iW-1][jW]);
            float rhov = 0.25*(W->rhov[iW][jW]+W->rhov[iW][jW-1]+W->rhov[iW-1][jW-1]+W->rhov[iW-1][jW]);
            float rhoE = 0.25*(W->rhoE[iW][jW]+W->rhoE[iW][jW-1]+W->rhoE[iW-1][jW-1]+W->rhoE[iW-1][jW]);
            float p = 0.25*(W->p[iW][jW]+W->p[iW][jW-1]+W->p[iW-1][jW-1]+W->p[iW-1][jW]);

            if (j==0 && i==0)
            {
                CellVolume = mesh->CellVolume[i][j];
            }
            else if (i==0 && j==mesh->nj-1)
            {
                CellVolume = mesh->CellVolume[i][j-1];
            }
            else if (i==mesh->ni-1 && j==0)
            {
                CellVolume = mesh->CellVolume[i-1][j];
            }
            else if(i==mesh->ni-1 && j==mesh->nj-1)
            {
                CellVolume = mesh->CellVolume[i-1][j-1];
            }
            else if (i==0)
            {
                CellVolume = 0.5*(mesh->CellVolume[i][j-1] + mesh->CellVolume[i][j]);
            }
            else if (i==mesh->ni-1)
            {
                CellVolume = 0.5*(mesh->CellVolume[i-1][j-1] + mesh->CellVolume[i-1][j]);
            }
            else if(j==0)
            {
                CellVolume = 0.5*(mesh->CellVolume[i][j] + mesh->CellVolume[i-1][j]);
            }
            else if (j==mesh->nj-1)
            {
                CellVolume = 0.5*(mesh->CellVolume[i][j-1] + mesh->CellVolume[i-1][j-1]);

            }
            else 
            {
                CellVolume = 0.25*(mesh->CellVolume[i][j]+mesh->CellVolume[i][j-1]+mesh->CellVolume[i-1][j-1]+mesh->CellVolume[i-1][j]);
            }
            fprintf(file, "%f %f %f %d %d %f %f %f %f %f\n", mesh->coordinates[i][j].x, mesh->coordinates[i][j].y, CellVolume, i, j, rhou, rhov, rho, rhoE, p);
        }
    }
    fclose(file);   
}