#ifndef MCLIB_H
#define MCLIB_H

#include "vector_operations.h"
#include "parameters.h"
#include "sample.h"

void AttemptParticleDisplacement (vec p[], int nP, int *accepted);
void AttemptParticleDeletion (vec p[], int *nP, int *accepted);
void AttemptParticleDeletionUS (vec p[], int *nP, int *accepted, int Nmin, int Nmax);
void AttemptParticleInsertion (vec p[], int *nP, int *accepted);
void AttemptParticleInsertionUS (vec p[], int *nP, int *accepted, int Nmin, int Nmax);
void AttemptParticleExchange (vec p[], int *nP, int *count_del, int *count_ins);
void AttemptParticleExchangeUS (vec p[], int *nP, int *count_del, int *count_ins, int Nmin, int Nmax);
void HistogramOfN (double H[], double R[], double HR[], double P[], FILE *prob);
void HistogramOfN_rec (double H[], double F[], double P[], FILE *prob);
void FCCInitialization(vec p[], int *nP);
void SCInitialization(vec p[], int *nP);
void BCCInitialization(vec p[], int *nP);
void InitRandomPosition (vec p[], int nP);
void InitConfiguration (FILE*point, vec p[], int *nP);
void SaveLastConfiguration (vec p[], int nP, int Nmin, int Nmax, int run);
void SaveLastSimulation (double H[], double F[], int run);
void UploadLastConfiguration (vec p[], int nP, int Nmin, int Nmax, int run);
void UploadLastSimulation (double H[], double F[], int run);


#endif
