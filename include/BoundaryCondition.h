#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "ConservedVariables.h"
#include "mesh.h"

// Enumeration for boundary condition types
typedef enum {
    WALL,
    FARFIELD,
    CONNECT
} BoundaryConditionType;

// Enumeration for boundary directions
typedef enum {
    I_DIRECTION,  
    J_DIRECTION  
} BoundaryDirection;

// Structure for boundary conditions
typedef struct {
    BoundaryConditionType type;       // Type of boundary condition
    BoundaryDirection direction;      // Direction of the boundary condition (along i or j information)
    int index_i;                      // Index in the i direction (set to -1 if not applicable)
    int index_j;                      // Index in the j direction (set to -1 if not applicable)
    // void (*apply)(ConservedVariables* W, MESH* mesh, int index_i, int index_j); // Function pointer to apply the boundary condition
} BoundaryCondition;



// Function to create boundary conditions
BoundaryCondition* createBoundaryCondition(BoundaryConditionType type, BoundaryDirection direction, int index_i, int index_j);
void applyWallBoundary(ConservedVariables* W, MESH* mesh, BoundaryCondition* bc_wall) ;


#endif // BOUNDARY_CONDITION_H
