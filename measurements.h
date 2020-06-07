#ifndef MEASURES_H
#define MEASURES_H
#include "vector_operations.h"
#include "cell_functions.h"
#include "physics.h"
#include "parameters.h"

void PrintAllPos (mol p[]);
void PrintPosFile(FILE* p, mol particles[], double t);
void PrintEnergyFile (FILE* p, double t, double kin_energy, double pot_energy);
void PrintMomFile (FILE* p, double t, vec sumv);
void PrintPressureFile (FILE* p, double t, double pressure);
double KinEnergy (mol particles[]);
double PotEnergy (mol particles[]);
void VelocityDistribution (FILE *p, FILE *max, mol particles[], double t_eq, double *res, int S);
double MaxwellBoltzmannDistribution (double v, double t_eq);
void EvaluateRDF (FILE*rdf, mol particles[], int sw, int* NRDF, double g[]);
void AndersenThermostat (mol particles[], int i, double T);
void Distances(double d[27], vec r1, vec r2);
void ModifyTemperature (double *T);
void EnergyAndPressureFromRDF (double g[], double *pressure, double *pot_energy, double t_eq);
void EvaluateTemperature (int sw, double pressure, double kin_energy, int *NT, double *t_eq, double *P);
void PrintVelocityField (mol particles[]);
void Density(mol particles[]);
void EvaluateInstantaneousRDF(FILE* rdf, mol particles[]);


#endif
