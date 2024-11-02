#ifndef WRITER_H
#define WRITER_H

#include "mesh.h"
#include "ConservedVariables.h"

void WriteTecplotFile(const char* filename, MESH* mesh, ConservedVariables* W);

#endif