#ifndef UTILS_H
#define UTILS_H


#define PI 3.14159265359

float L2norm(float* P1,float* P2);
float* calculateNormalVector(float* P1, float* P2, float ds);
double getCellVolume(float* a, float*b);
void display_image(const char *image_path);
#endif
