#ifndef SAMPLE_H
#define SAMPLE_H

#include "vector_operations.h"
// #include "cell_functions.h"
// #include "physics.h"
#include "parameters.h"

double Energy (vec p[], int nP, vec pos, int k);
void PrintAllPosFile(FILE* point, vec particles[], int nP, int step);
void PrintAllPos(vec particles[], int nP, int step);
void Info (int nPin, int nPfin, double time, double P[]);
void SampleNumberOfParticles(FILE *num, int nP, int i);
double TotEnergy (vec particles[], int nP);
void PrintEnergyFile (FILE* en, int i, vec p[], int nP);
void FreeEnergy (double P[], FILE *en);
double Pressure (double P0);
double AverageDensity (double P[]);
void HistogramReweightingZ (double new_zeta, double P[], double newP[]);


#endif
