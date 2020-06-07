#ifndef VECTOR_H
#define VECTOR_H

#include "parameters.h"

typedef struct {
  double x, y, z;
} vec;

typedef struct {
  vec r, v;
} mol;

void PrintD (double x, char*);
void PrintI (int x, char*);
void PrintParticle (mol a);
void PrintIntArray (char* name, int v[], int l);
void PrintRealArray (char* name, double v[], int l);
void PrintVecAndName (vec, char*);
void PrintStruct (mol);
void PrintVec (vec);

vec PeriodicImage (vec r);
vec NearestIntVec (vec v);
vec SeparationPBC (vec r_i, vec r_j);
double Distance (vec r_i, vec r_j);

double RandReal(double a, double b);
int RandInt(int a, int b);
double GaussianNumber (double sigma, double mean);

double Cube (double a);

vec SetVec (double, double, double);
vec DivVec (vec , double);
vec SumVec (vec, vec);
vec SubVec (vec , vec );
double ScalarProduct (vec a, vec b);
double SquaredNorm (vec);
vec MultiplyVec (vec, double);
vec SetEqualVec (vec a);
mol SetEqualMol (mol a);
double Mod (double a, double b);
void SetZeroRealArray(double v[], int length);
void SetOneRealArray (double v[], int length);
void SetOneIntArray(int v[], int length);
void SetZeroIntArray(int v[], int length);
void SetZeroArrayVec (vec v[], int length);
vec SetConstantVector (double c);

void Space();

#endif
