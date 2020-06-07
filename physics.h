#ifndef PHYSICS_H
#define PHYSICS_H
#include "vector_operations.h"
#include "cell_functions.h"
#include "parameters.h"
#include "measurements.h"

void Info (mol particles[]);
void InitPosition (mol particles[]);
void InitVel (mol particles[]);
void Setup (mol particles[], double *kin_energy, double *pot_energy);
double LJPotential(double r);
double DerivativeLJ(double r);
vec Force (vec r_i, vec r_j, double *pot_energy_couple, double *p);
void Neighbours (int neigh[1+nCells][14], vec cell);
void BuildNeighbours (int neigh[1+nCells][14]);
void MoveParticle (mol particles[], vec displ[], vec f, int i, double T, double t, double EQ_TIME);
void FirstLeapfrogStep (mol particles[], int cellnumber[], int neigh[1+nCells][14], int head[1+nCells], int after[nParticles+1], int before[nParticles+1], double *pressure);
void TakeStep (int cellnumber[], mol particles[], int head[], int after[], int neigh[1+nCells][14], double *pot_energy, double *kin_energy, vec* q, double *pressure, double *msd, vec displ[], double T, double t, double EQ_TIME);
void AlternativeInitPosition (mol particles[]);


#endif
