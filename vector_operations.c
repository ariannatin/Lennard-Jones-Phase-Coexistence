#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector_operations.h"

void PrintVec (vec a) {
  printf("(%lf, %lf, %lf)\n", a.x, a.y, a.z);
}

void PrintVecAndName (vec a, char* s) {
  printf("%s (%lf, %lf, %lf)\n", s, a.x, a.y, a.z);
}

void PrintIntArray(char* name, int v[], int l) {
  int i;
  printf("%s  ", name);
  for(i=0; i<l; i++) {
    printf("%d ", v[i]);
  }
  Space();
}

void PrintRealArray (char* name, double v[], int l) {
  int i;
  printf("%s  ", name);
  for(i=0; i<l; i++) {
    printf("%lf ", v[i]);
  }
  Space();
}

void PrintParticle (mol a) {
  printf("Position: ");
  PrintVec(a.r);
  printf("Velocity: ");
  PrintVec(a.v);
}

vec PeriodicImage (vec r) {
  r.x += 0.5*L;
  r.x = Mod(r.x, L);
  r.x -= 0.5*L;

  r.y += 0.5*L;
  r.y = Mod(r.y, L);
  r.y -= 0.5*L;

  r.z += 0.5*L;
  r.z = Mod(r.z, L);
  r.z -= 0.5*L;

  return r;
}

vec NearestIntVec(vec v) {
  vec r;
  r.x = round(v.x);
  r.y = round(v.y);
  r.z = round(v.z);
  return r;
}

vec SeparationPBC (vec r_i, vec r_j) {
  vec r_ij = SubVec(r_i, r_j), R;
  R = SubVec(r_ij, MultiplyVec(NearestIntVec(MultiplyVec(r_ij, Li)), L));
  return R;
}

double Distance (vec r_i, vec r_j) {
  return  sqrt(SquaredNorm(SeparationPBC(r_i, r_j)));
}

double SquaredNorm (vec r) {
  return (r.x * r.x + r.y * r.y + r.z * r.z);
}

double Cube (double a) {
  return a*a*a;
}

vec SumVec (vec a, vec b) {
  return SetVec (a.x + b.x, a.y + b.y, a.z + b.z);
}

vec SubVec (vec a, vec b) {
  return SetVec (a.x - b.x, a.y - b.y, a.z - b.z);
}

vec MultiplyVec (vec a, double t) {
  return SetVec (a.x * t, a.y * t, a.z * t);
}

vec DivVec (vec a, double t) {
  return SetVec(a.x / t, a.y / t, a.z / t);
}

double ScalarProduct (vec a, vec b) {
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

vec SetVec (double a, double b, double c) {
  vec v;
  v.x = a;
  v.y = b;
  v.z = c;
  return v;
}

vec SetEqualVec (vec a) {
  return SetVec (a.x, a.y, a.z);
}

mol SetEqualMol (mol a) {
  mol m;
  m.r = SetEqualVec (a.r);
  m.v = SetEqualVec (a.v);
  return m;
}

double RandReal (double a, double b) {
  return drand48()*(b-a) + a;
}

int RandInt(int a, int b){
  return rand()/(RAND_MAX +1.) * (b+1-a) +a ;
}

double Mod(double a, double b) {
  double r = fmod(a, b);
  return r < 0 ? r + b : r;
}

double GaussianNumber(double sigma, double mean) {
  double r = 2., v1, v2, l;
  do {
    v1 = 2. * RandReal(0., 1.) - 1.;
    v2 = 2. * RandReal(0., 1.) - 1.;
    r = v1 * v1 + v2 * v2;
  } while(r>=1.);
  l = v1 * sqrt(-2.*log(r)/r);
  l = mean + sigma * l;
  return l;
}

void SetZeroRealArray(double v[], int length) {
  int i;
  for(i=0; i<length; i++) {
    v[i] = 0.;
  }
}

void SetOneRealArray (double v[], int length) {
  int i;
  for(i=0; i<length; i++) {
    v[i] = 1.;
  }
}

void SetZeroIntArray(int v[], int length) {
  int i;
  for(i=0; i<length; i++) {
    v[i] = 0;
  }
}

void SetOneIntArray(int v[], int length) {
  int i;
  for(i=0; i<length; i++) {
    v[i] = 1;
  }
}

void SetZeroArrayVec (vec v[], int length) {
  int i;
  for(i=0; i<length; i++) {
    v[i] = SetVec(0., 0., 0.);
  }
}

vec SetConstantVector (double c) {
  return SetVec(c, c, c);
}

void PrintD (double x, char* name) {
  printf("%s: %.20lf\n", name, x);
}

void PrintI (int x, char* name) {
  printf("%s: %d\n", name, x);
}

void Space() {
  printf("\n");
}
