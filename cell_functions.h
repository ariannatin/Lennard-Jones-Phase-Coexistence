#ifndef CELL_H
#define CELL_H

#include "vector_operations.h"
#include "parameters.h"

int NumberOfCellFromPosition (vec r);
int NumberOfCell (vec v);
void CreateCellList (int cellnumber[], int head[], int after[], int before[], mol particles[]);
void UpdateCellParticle (int head[], int after[], int before[], int i, int old_icell, int new_icell);
void UpdateCellList (int head[], int after[], int before[], int cellnumber[], mol particles[]);

#endif
