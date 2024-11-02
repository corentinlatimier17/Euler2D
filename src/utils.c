#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265359

float L2norm(float* P1,float* P2)
{
    return sqrt(pow(P1[0]-P2[0],2)+pow(P1[1]-P2[1],2));
}

// Function to calculate and return the tangent vector as a dynamically allocated array
float* calculateNormalVector(float* P1, float* P2, float ds) {
    float dx = P2[0] - P1[0];
    float dy = P2[1] - P1[1];

    // Allocate memory for the tangent vector
    float* tangent = (float*)malloc(2 * sizeof(float));
    if (tangent == NULL) {
        fprintf(stderr, "Memory allocation failed for tangent vector\n");
        exit(1);  // Handle memory allocation failure
    }

    // Calculate unit tangent vector
    if (ds != 0) {
        tangent[0] = dx / ds;
        tangent[1] = dy / ds;
    } else {
        tangent[0] = 0;
        tangent[1] = 0;
    }

    // Allocate memory for the normal vector
    float* normal = (float*)malloc(2 * sizeof(float));
    if (normal == NULL) {
        fprintf(stderr, "Memory allocation failed for normal vector\n");
        exit(1);  // Handle memory allocation failure
    }

    // Apply 90° rotation (anticlockwise direction)
    normal[0] = -tangent[1];
    normal[1] = tangent[0]; 

    free(tangent);
    return normal;
}

double getCellVolume(float* a, float* b)
{
    return 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
}

void display_image(const char *image_path) {
    system("figlet  -w 100 Euler   2D   CFD  solver");
    char command[256];  // Assurez-vous que la taille est suffisante pour le chemin
    snprintf(command, sizeof(command), "jp2a --width=90 --height=35 --invert %s", image_path);
    system(command);
}