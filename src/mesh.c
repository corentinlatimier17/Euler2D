#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "mesh.h"

void ReadCoordinates(MESH* mesh)
{
    FILE* file = fopen(mesh->filepath, "r");
    if (file == NULL) {
        perror("Error opening mesh file");
        exit(EXIT_FAILURE);
    }

    int num_blocks = 0;
    int ni = 0;
    int nj = 0;

    // Read the first line (number of blocks)
    if (fscanf(file, "%d", &num_blocks) != 1) {
        fprintf(stderr, "Error reading number of blocks.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Read second line (ni, nj)
    if (fscanf(file, "%d %d", &ni, &nj) != 2) {
        fprintf(stderr, "Error reading ni and nj.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Read x coordinates first (order: i,j)
    for (int i = 0; i < mesh->ni; i++)
    {
        for (int j = 0; j < mesh->nj; j++)
        {
            if (fscanf(file, "%f", &mesh->coordinates[i][j].x) != 1) {
                fprintf(stderr, "Error reading x coordinate at (%d, %d).\n", i, j);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }
    // Read y coordinates then (order: i,j)
    for (int i = 0; i < mesh->ni; i++)
    {
        for (int j = 0; j < mesh->nj; j++)
        {
            if (fscanf(file, "%f", &mesh->coordinates[i][j].y) != 1) {
                fprintf(stderr, "Error reading y coordinate at (%d, %d).\n", i, j);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    fclose(file);
}

MESH* LoadMESH(char* filepath)
{
    MESH* mesh = (MESH*)malloc(sizeof(MESH));

    if (mesh == NULL) {
        perror("Failed to allocate memory for mesh structure");
        return NULL; // Handle memory allocation failure
    }

    // Set the filepath attribute to the object
    mesh->filepath = filepath;
    
    // Read plot3d file header (nblocks, ni, nj)
    ReadHeaderPlot3D(mesh);

    // Allocate memory for the coordinates array based on ni and nj
    mesh->coordinates = (VerticeCoordinate**)malloc(mesh->ni * sizeof(VerticeCoordinate*));
    if (mesh->coordinates == NULL) {
        perror("Failed to allocate memory for mesh->coordinates matrix");
        free(mesh); // Free the allocated mesh structure before returning
        return NULL; // Handle memory allocation failure
    }

    for (int i = 0; i < mesh->ni; i++)
    {
        mesh->coordinates[i] = (VerticeCoordinate*)malloc(mesh->nj * sizeof(VerticeCoordinate));
        if (mesh->coordinates[i] == NULL) {
            perror("Failed to allocate memory for one row in mesh->coordinates");
            // Free previously allocated rows
            for (int k = 0; k < i; k++) {
                free(mesh->coordinates[k]);
            }
            free(mesh->coordinates);
            free(mesh);
            return NULL; // Handle memory allocation failure
        }
    }

    ReadCoordinates(mesh); // fill mesh.coordinates with x and y values in a structured way


    // Allocate memory for the vertical faces storage
    mesh->VerticalFace = (Face**)malloc(sizeof(Face*)*(mesh->ni-1));
    if (mesh->VerticalFace == NULL) 
        {
        perror("Failed to allocate memory for the matrix of vertical faces");
        exit(EXIT_FAILURE);
        }
    for (int i = 0; i < (mesh->ni-1); i++) 
    {
        mesh->VerticalFace[i] = (Face*)malloc(mesh->nj * sizeof(Face));
        if (mesh->VerticalFace[i] == NULL) {
            perror("Failed to allocate memory for one row of the matrix of vertical faces");
            exit(EXIT_FAILURE);
        }
    }
    
    // Allocate memory for the horizontal faces storage
    mesh->HorizontalFace = (Face**)malloc(sizeof(Face*)*(mesh->ni));
    if (mesh->HorizontalFace == NULL) 
        {
        perror("Failed to allocate memory for the matrix of horizontal faces");
        exit(EXIT_FAILURE);
        }
    for (int i = 0; i < (mesh->ni); i++) 
    {
        mesh->HorizontalFace[i] = (Face*)malloc((mesh->nj-1) * sizeof(Face));
        if (mesh->HorizontalFace[i] == NULL) 
        {
            perror("Failed to allocate memory for one row of the matrix of horizontal faces");
            exit(EXIT_FAILURE);
        }
    }

    // Initialization of faces
    InitializeFaces(mesh);

    mesh->CellVolume = (double**)malloc(sizeof(double*) * (mesh->ni-1));
    if (mesh->CellVolume == NULL) 
    {
        perror("Failed to allocate memory for mesh->CellVolume");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < (mesh->ni-1); i++) 
    {
        mesh->CellVolume[i] = (double*)malloc(sizeof(double) * (mesh->nj-1));
        if (mesh->CellVolume[i] == NULL) 
        {
            perror("Failed to allocate memory for a row in mesh->CellVolume");
        
            // Free previously allocated rows in case of failure
            for (int k = 0; k < i; k++) 
            {
                free(mesh->CellVolume[k]);
            }
            free(mesh->CellVolume);
            exit(EXIT_FAILURE);
        }
    }

    for (int i=0; i<mesh->ni-1; i++)
    {
        for (int j=0; j<mesh->nj-1; j++)
        {
            // Schematic of the points

            // P1 ******** P2
            //    *      *   
            //    *      *
            // P3 ******** P4   

            float* P1 = getIJtoXY(mesh, i,j);
            float* P2 = getIJtoXY(mesh, i, j+1);
            float* P3 = getIJtoXY(mesh, i+1,j);
            float* P4 = getIJtoXY(mesh, i+1, j+1);

            float* a = (float*)malloc(sizeof(float)*2);
            a[0] = P4[0]-P1[0];
            a[1] = P4[1]-P1[1];
            float* b = (float*)malloc(sizeof(float)*2);
            b[0] = P3[0]-P2[0];
            b[1] = P3[1]-P2[1];
            free(P1);
            free(P2);
            free(P3);
            free(P4);
            mesh->CellVolume[i][j] = getCellVolume(a, b);
            free(a);
            free(b);
        }
    }
    return mesh; // Return the pointer to the allocated mesh
}

void InitializeFaces(MESH* mesh) {
    for (int i = 0; i < mesh->ni - 1; i++) { // Vertical faces
        for (int j = 0; j < mesh->nj; j++) {
            float* P1 = getIJtoXY(mesh, i, j);
            float* P2 = getIJtoXY(mesh, i + 1, j);
            mesh->VerticalFace[i][j].ds = L2norm(P1, P2);
            float* normal = calculateNormalVector(P1, P2, mesh->VerticalFace[i][j].ds);
            mesh->VerticalFace[i][j].normal.nx = normal[0];
            mesh->VerticalFace[i][j].normal.ny = normal[1];
            free(normal);
            free(P1);
            free(P2);
        }
    }

    for (int i = 0; i < mesh->ni; i++) { // Horizontal faces
        for (int j = 0; j < mesh->nj - 1; j++) {
            float* P1 = getIJtoXY(mesh, i, j);
            float* P2 = getIJtoXY(mesh, i, j + 1);
            mesh->HorizontalFace[i][j].ds = L2norm(P1, P2);
            float* normal = calculateNormalVector(P1, P2, mesh->HorizontalFace[i][j].ds);
            mesh->HorizontalFace[i][j].normal.nx = normal[0];
            mesh->HorizontalFace[i][j].normal.ny = normal[1];
            free(normal);
            free(P1);
            free(P2);
        }
    }
}

void ReadHeaderPlot3D(MESH* mesh)
{
    FILE* file = fopen(mesh->filepath, "r");
    if (file == NULL) {
        perror("Error opening mesh file");
        exit(EXIT_FAILURE);
    }

    // Read the first line (number of blocks)
    fscanf(file, "%d", &mesh->num_blocks);
    // Read second line (ni, nj)
    fscanf(file, "%d %d", &mesh->ni, &mesh->nj);

    mesh->nVFaces = (mesh->ni-1)*mesh->nj; // number of vertical faces
    mesh->nHFaces = (mesh->nj-1)*mesh->ni; // number of horizontal faces

    fclose(file);
}

float* getIJtoXY(MESH* mesh, int i, int j) // get x,y coordinates for given (i,j) indexes 
{
    if (i<0 || i>=mesh->ni || j<0 || j>=mesh->nj)
    {
        perror("Indexes i or j provided are invalid");
        exit(EXIT_FAILURE);
    }
    else
    {
        // Allocate a 2-element array
        float* coordinates = (float*)malloc(2 * sizeof(float));
        if (coordinates == NULL) {
            perror("Failed to allocate memory for coordinates");
            return NULL;
        }

        // Assign x and y to the array
        coordinates[0] = mesh->coordinates[i][j].x;
        coordinates[1] = mesh->coordinates[i][j].y;

        return coordinates;
    }
}


Normal2D getNormal(int i1, int j1, int i2, int j2, MESH* mesh)
{
    if (i1<0 || i2<0 || j1<0 || j2<0 || i1>=mesh->ni || i2>=mesh->ni || j1>=mesh->nj || j2>= mesh->nj)
    {
        perror("getNormal : invalid indexes (i,j)");
        exit(EXIT_FAILURE);
    }
    if ((abs(i1-i2)==1 && abs(j1-j2)==0) || (abs(j1-j2)==1 && abs(i1-i2)==0))
    { 
        if (abs(i1-i2)==1) // vertical face
        {
            int i = fmin(i1,i2);
            return mesh->VerticalFace[i][j1].normal;
        }
        if (abs(j1-j2)==1) // horizontal face
        {
            int j = fmin(j1,j2);
            return mesh->HorizontalFace[i1][j].normal;
        }
    }
    else
    {
        perror("getNormal : points are not connected with a face");
        exit(EXIT_FAILURE);
    }
}


float getDS(int i1, int j1, int i2, int j2, MESH* mesh)
{
    if (i1<0 || i2<0 || j1<0 || j2<0 || i1>=mesh->ni || i2>=mesh->ni || j1>=mesh->nj || j2>= mesh->nj)
    {
        perror("getDS : invalid indexes (i,j)");
        exit(EXIT_FAILURE);
    }
    if ((abs(i1-i2)==1 && abs(j1-j2)==0) || (abs(j1-j2)==1 && abs(i1-i2)==0))
    { 
        if (abs(i1-i2)==1) // vertical face
        {
            int i = fmin(i1,i2);
            return mesh->VerticalFace[i][j1].ds;
        }
        if (abs(j1-j2)==1) // horizontal face
        {
            int j = fmin(j1,j2);
            return mesh->HorizontalFace[i1][j].ds;
        }
    }
    else
    {
        perror("getDS : points are not connected with a face");
        exit(EXIT_FAILURE);
    }
}

// Function to free the allocated memory for a mesh
void FreeMESH(MESH* mesh) {
    if (mesh != NULL) {
        // Free the coordinates matrix
        if (mesh->coordinates != NULL) {
            for (int i = 0; i < mesh->ni; i++) {
                free(mesh->coordinates[i]); // Free each row
            }
            free(mesh->coordinates); // Free the array of row pointers
        }

        // Free the vertical normals matrix
        if (mesh->VerticalFace != NULL) {
            for (int i = 0; i < mesh->ni - 1; i++) {
                free(mesh->VerticalFace[i]); // Free each row
            }
            free(mesh->VerticalFace); // Free the array of row pointers
        }

        // Free the horizontal normals matrix
        if (mesh->HorizontalFace != NULL) {
            for (int i = 0; i < mesh->ni; i++) {
                free(mesh->HorizontalFace[i]); // Free each row
            }
            free(mesh->HorizontalFace); // Free the array of row pointers
        }

        // Free the mesh structure itself
        free(mesh);
    }
}


void PrintYCoordinates(MESH* mesh) {
    printf("Y Coordinates:\n");
    for (int i = 0; i < mesh->ni; i++) {
        for (int j = 0; j < mesh->nj; j++) {
            printf("%f ", mesh->coordinates[i][j].y);
            if (j== mesh->nj-1)
            {
                printf("\n");
            }
        }
    }
}

void PrintXCoordinates(MESH* mesh) {
    printf("X Coordinates:\n");
    for (int i = 0; i < mesh->ni; i++) {
        for (int j = 0; j < mesh->nj; j++) {
            printf("%f ", mesh->coordinates[i][j].x);
            if (j== mesh->nj-1)
            {
                printf("\n");
            }
        }
    }
}

/// Verification of the metric 
float*** allocateResidual(int ni, int nj) {
    if (ni <= 1 || nj <= 1) {
        fprintf(stderr, "Invalid mesh dimensions for allocation.\n");
        return NULL;
    }
    float*** residual = (float***)calloc(ni - 1, sizeof(float**));
    if (!residual) return NULL; // Check if allocation failed
    
    for (int i = 0; i < ni - 1; i++) {
        residual[i] = (float**)calloc(nj - 1, sizeof(float*));
        if (!residual[i]) {
            // Free previously allocated rows before returning
            for (int k = 0; k < i; k++) free(residual[k]);
            free(residual);
            return NULL;
        }
        for (int j = 0; j < nj - 1; j++) {
            residual[i][j] = (float*)calloc(2, sizeof(float)); // Allocate and initialize to 0
            if (!residual[i][j]) {
                // Free previously allocated rows before returning
                for (int k = 0; k <= i; k++) {
                    for (int l = 0; l < (k == i ? j : nj - 1); l++) {
                        free(residual[k][l]);
                    }
                    free(residual[k]);
                }
                free(residual);
                return NULL;
            }
        }
    }
    return residual;
}

void freeResidual(float*** residual, int ni, int nj) {
    for (int i = 0; i < ni - 1; i++) {
        for (int j = 0; j < nj - 1; j++) {
            free(residual[i][j]);
        }
        free(residual[i]);
    }
    free(residual);
}

void verifyMetric(MESH* mesh) {
    float*** residual = allocateResidual(mesh->ni, mesh->nj);
    if (!residual) {
        fprintf(stderr, "Memory allocation for residual failed.\n");
        exit(EXIT_FAILURE);
    }

    // 1. Loop on horizontal faces
    for (int i = 0; i < mesh->ni; i++) {
        for (int j = 0; j < mesh->nj - 1; j++) {
            float s = getDS(i, j, i, j + 1, mesh);
            Normal2D normal = getNormal(i, j, i, j + 1, mesh);

            if (i == 0) {
                residual[i][j][0] += s * normal.nx;
                residual[i][j][1] += s * normal.ny;
            } else if (i == mesh->ni - 1) {
                residual[i - 1][j][0] -= s * normal.nx; // Adjusted index for proper access
                residual[i - 1][j][1] -= s * normal.ny;
            } else {
                residual[i][j][0] += s * normal.nx;
                residual[i][j][1] += s * normal.ny;
                residual[i - 1][j][0] -= s * normal.nx;
                residual[i - 1][j][1] -= s * normal.ny;
            }
        }
    }

    // 2. Loop on vertical faces
    for (int i = 0; i < mesh->ni - 1; i++) {
        for (int j = 0; j < mesh->nj; j++) {
            float s = getDS(i, j, i + 1, j, mesh);
            Normal2D normal = getNormal(i, j, i + 1, j, mesh);

            if (j == 0) {
                residual[i][j][0] -= s * normal.nx;
                residual[i][j][1] -= s * normal.ny;
            } else if (j == mesh->nj - 1) {
                residual[i][j - 1][0] += s * normal.nx; // Adjusted index for proper access
                residual[i][j - 1][1] += s * normal.ny;
            } else {
                residual[i][j][0] -= s * normal.nx;
                residual[i][j][1] -= s * normal.ny;
                residual[i][j - 1][0] += s * normal.nx;
                residual[i][j - 1][1] += s * normal.ny;
            }
        }
    }

    // Print the residual values
    for (int i = 0; i < mesh->ni - 1; i++) {
        for (int j = 0; j < mesh->nj - 1; j++) {
            printf("(%f, %f) ", residual[i][j][0], residual[i][j][1]);
        }
        printf("\n"); // Newline for the next row
    }

    // Free residual memory at the end
    freeResidual(residual, mesh->ni, mesh->nj);
}
