#ifndef GRID_H
#define GRID_H

typedef struct{
    float x;
    float y;
}VerticeCoordinate;


typedef struct {
    float nx;
    float ny;
} Normal2D;

typedef struct
{
    Normal2D normal;
    float ds;
}Face;


typedef struct{
    int ni;
    int nj;
    int nVFaces;
    int nHFaces;
    int num_blocks;
    VerticeCoordinate** coordinates;
    double** CellVolume;
    Face** VerticalFace;
    Face** HorizontalFace;
    char* filepath;
}MESH;


void ReadCoordinates(MESH* mesh);
MESH* LoadMESH(char* filepath);
void InitializeFaces(MESH* mesh);
void ReadHeaderPlot3D(MESH* mesh);
float* getIJtoXY(MESH* mesh, int i, int j);
Normal2D getNormal(int i1, int j1, int i2, int j2, MESH* mesh);
float getDS(int i1, int j1, int i2, int j2, MESH* mesh);
void FreeMESH(MESH* mesh);
void PrintYCoordinates(MESH* mesh);
void PrintXCoordinates(MESH* mesh);

// for verification of the implementation
float*** allocateResidual(int ni, int nj);
void freeResidual(float*** residual, int ni, int nj);
void verifyMetric(MESH* mesh);

#endif
