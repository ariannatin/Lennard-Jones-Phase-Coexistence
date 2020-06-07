#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sample.h"

double Energy (vec p[], int nP, vec pos, int k) {
  int i;
  double u = 0., r2, r2i, r6i;
  for(i=0; i<nP; i++) {
    if(i!=k) {
      r2 = SquaredNorm(SeparationPBC(p[i], pos));
      if(r2<r2_cut) {
        r2i = 1./r2;
        r6i = Cube(r2i);
        u += 4 * r6i * (r6i - 1);
        // u += 0.;
      }
    }
  }
  return u;
}

double TotEnergy (vec particles[], int nP) {
  int i, j;
  double u = 0., r2, r2i, r6i;
  for(i=0; i<nP; i++) {
    for(j=(i+1); j<nP; j++) {
      r2 = SquaredNorm(SeparationPBC(particles[i], particles[j]));
      if(r2<r2_cut) {
        r2i = 1./r2;
        r6i = Cube(r2i);
        u += 4 * r6i * (r6i - 1);
      }
    }
  }
  return u/(double)nP;
}

void PrintEnergyFile (FILE* en, int i, vec p[], int nP) {
  fprintf(en, "%d %lf\n", i, TotEnergy(p, nP));
}


void PrintAllPosFile(FILE* point, vec particles[], int nP, int step) {
  int i;
  for(i=0; i<nP; i++) {
    // fprintf(point, "%d %.20lf %.20lf %.20lf\n", step, particles[i].x, particles[i].y, particles[i].z);
    fprintf(point, "%.20lf %.20lf %.20lf\n", particles[i].x, particles[i].y, particles[i].z);
  }
  fprintf(point, "\n\n");
}

void PrintAllPos(vec particles[], int nP, int step) {
  int i;
  for(i=0; i<nP; i++) {
    // printf("%d %.20lf %.20lf %.20lf\n", step, particles[i].x, particles[i].y, particles[i].z);
    printf("particella %d %.20lf %.20lf %.20lf\n", i, particles[i].x, particles[i].y, particles[i].z);
  }
}

void Info (int nPin, int nPfin, double time, double P[]) {
  printf("\n\n*******************************************************\n");
  printf("      This is a simulation of a Lennard-Jones fluid\n             with the following parameters:");
  printf("\n=======================================================\n");
  printf("        Initial particles: %d          \n", nPin);
  printf("        Final particles: %d          \n", nPfin);
  printf("        Initial Density: %lf                 \n", (double)nPin/(double)V);
  printf("        Final Density: %lf                 \n", (double)nPfin/(double)V);
  printf("        Cube of %.4lf x %.4lf x %.4lf                 \n", L, L, L);
  Space();
  // printf("      Number of cells per side: %d            \n", W);
  printf("        Temperature: %.2lf                \n", T);
  printf("        Fugacity: %.5lf                \n", zeta);
  Space();
  printf("        Monte Carlo steps: %d                \n", nSteps);
  printf("        Exchange attempts per step: %d                \n", nEx);
  printf("        Displacement attempts per step: %d                \n", nDispl);
  Space();
  printf("        Average Density: %lf                 \n", AverageDensity(P));
  printf("        Pressure: %lf                      \n", Pressure(P[0]));
  Space();
  printf("        Execution time:\n              %lf s\n              %lf min\n              %lf h\n", time, time/60., time/3600.);
  printf("=======================================================\n");
}

void SampleNumberOfParticles(FILE *num, int nP, int i) {
  fprintf(num, "%d %d\n", i, nP);
}

void HistogramReweightingZ (double new_zeta, double P[], double newP[]) {
  int n;
  double sum = 0.;
  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    newP[n] = pow((new_zeta/zeta), n)*P[n];
    sum += newP[n];
  }
  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    newP[n] /= sum;
  }
}

double Pressure (double P0) {
  if (P0 == 0) {
    return 0;
  } else {
    return -T/V*log(P0);
  }
}

double AverageDensity (double P[]) {
  int n;
  double sum = 0.;
  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    sum += n*P[n]/V;
  }
  return sum;
}

void FreeEnergy (double P[], FILE *en) {
  int n;
  double A[N_SAMPLE_MAX+1];

  for(n=0; n<(N_SAMPLE_MAX+1); n++) {
    if(P[0] >= 1e-100) {
      A[n] = n*T*log(zeta) - T*log(P[n]/P[0]);
      fprintf(en, "%lf\n", A[n]);
    }
  }
}
