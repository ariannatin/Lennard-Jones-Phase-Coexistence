#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "measurements.h"

void PrintAllPos (mol p[]) {
  int i;
  for (i=0; i<nParticles; i++) {
    PrintParticle(p[i]);
    printf("\n");
  }
}

void PrintPosFile(FILE* p, mol particles[], double t) {
  int i;
  for(i=0; i<nParticles; i++)  {
    fprintf(p, "%lf %.20lf %.20lf %.20lf\n", t, particles[i].r.x, particles[i].r.y, particles[i].r.z);
  }
  // fprintf(p, "%lf %.20lf %.20lf %.20lf\n", t, particles[1].r.x, particles[1].r.y, particles[1].r.z);
  fprintf(p, "\n\n");
}

void PrintEnergyFile (FILE* p, double t, double kin_energy, double pot_energy) {
  double T = 2.*kin_energy/DIM;
  fprintf(p, "%.10lf %.20lf %.20lf %.20lf %.10lf\n", t, kin_energy, pot_energy, pot_energy+kin_energy, T);
}

void PrintMomFile (FILE* p, double t, vec sumv) {
  fprintf(p, "%lf %.20lf %.20lf %.20lf\n", t, sumv.x, sumv.y, sumv.z);
}

void PrintPressureFile (FILE* p, double t, double pressure) {
  fprintf(p, "%lf %.20lf\n", t, pressure);
}

double KinEnergy (mol particles[]) {
  int i;
  double k = 0;
  for(i=0; i<nParticles; i++) {
    k += SquaredNorm(particles[i].v);
  }
  return k/(2.*nParticles);
}

double PotEnergy (mol particles[]) {
  int i, j;
  double u = 0., r2, r2i, r6i;
  for(i=0; i<nParticles; i++) {
    for(j=(i+1); j<nParticles; j++) {
      r2 = SquaredNorm(SeparationPBC(particles[i].r, particles[j].r));
      if(r2<r2_cut) {
        r2i = 1./r2;
        r6i = Cube(r2i);
        u += 4 * r6i * (r6i - 1) + 4*e_cut;
      }
    }
  }
  return u/(double)nParticles;
}

void VelocityDistribution (FILE *p, FILE *max,  mol particles[], double t_eq, double *res, int S) {
  int i;
  double xmax = -1.;
  double vv[nParticles];

  // printf("Starting computing velocity distribution...\n");

  double meanV = 0.;

  for(i=0; i<nParticles; i++){
    vv[i] = sqrt(SquaredNorm(particles[i].v));
    meanV += vv[i];
    if(vv[i]>=xmax) xmax = vv[i];
  }

  meanV /= nParticles;

  int n;
  int nBins = xmax/binsize + 1;

  double *h = (double*) calloc(nBins, sizeof(double));

  if(h == NULL) {
    printf("Allocation failed\n");
    exit(0);
  }

  for(n=0; n<nParticles; n++) {
    i =  vv[n]/binsize;
    h[i] = h[i] + 1.;
  }

  *res = 0.;

  double xx, mb;
  for(n=0; n<nBins; n++) {
    h[n] /= (binsize*nParticles);
    xx = (n+0.5)*binsize;
    if(S==1) fprintf(p, "%lf %lf\n", xx, h[n]);
    printf("%lf %lf\n", xx, h[n]);
    mb = MaxwellBoltzmannDistribution(xx, t_eq);
    if(S==1) fprintf(max, "%lf %lf\n", xx, mb);
    *res += (h[n] - mb)*(h[n] - mb)/mb;
  }

  // *res /= (double)nBins;
  // printf("Done!\n");
  free(h);

}

double MaxwellBoltzmannDistribution (double v, double t_eq) {
  double a = 4*M_PI*v*v*exp(-v*v/(2*t_eq))*pow(1./(2*t_eq*M_PI), 1.5);

  return a;
}

void EvaluateRDF (FILE*rdf, mol particles[], int sw, int* NRDF, double g[]) {
  #ifdef RDF

  double r, vb, nid;
  int i, j, n, m;
  double d[27];

  switch (sw) {
    case 0:

      if((*NRDF)==0) printf("Starting computing Radial Distribution Function...\n");

      (*NRDF)++;
      for(i=0; i<nParticles; i++) {
        for(j=i+1; j<nParticles; j++) {
          // r = sqrt(SquaredNorm(SubVec(particles[i].r, particles[j].r)));
          Distances(d, particles[i].r, particles[j].r);
          for(m=0; m<27; m++) {
            if(d[m]<L) {   // forse c'è un errore?
              n = d[m]/binsize_RDF;
              g[n] += 2;
            }
          }
        }
      }
      break;
    case 1:

      for(i=0; i<nBins_RDF; i++) {
        r = binsize_RDF*(i + 0.5);
        vb = (Cube((double)(i+1)) - Cube(i))*Cube(binsize_RDF);
        nid = 4./3*M_PI*vb*rho;
        g[i] = g[i]/((*NRDF)*nParticles*nid);
        fprintf(rdf, "%lf %.10lf\n", r, g[i]);
      }
      break;
  }

  #endif
}

void EvaluateInstantaneousRDF(FILE* rdf, mol particles[]) {
  double r, vb, nid;
  int i, j, n, m;
  double d[27];
  double g[nBins_RDF];
  SetZeroArray(g, nBins_RDF);

  for(i=0; i<nParticles; i++) {
    for(j=i+1; j<nParticles; j++) {
      Distances(d, particles[i].r, particles[j].r);
      for(m=0; m<27; m++) {
        if(d[m]<L) {
          n = d[m]/binsize_RDF;
          g[n] += 2;
        }
      }
    }
  }

  for(i=0; i<nBins_RDF; i++) {
    r = binsize_RDF*(i + 0.5);
    vb = (Cube((double)(i+1)) - Cube(i))*Cube(binsize_RDF);
    nid = 4./3*M_PI*vb*rho;
    g[i] = g[i]/(nParticles*nid);
    fprintf(rdf, "%lf %.10lf\n", r, g[i]);
  }
}

void AndersenThermostat (mol particles[], int i, double T) {
  if(RandReal(0., 1.) < FREQ*dt) {
    particles[i].v.x = GaussianNumber(sqrt(T), 0.);
    particles[i].v.y = GaussianNumber(sqrt(T), 0.);
    particles[i].v.z = GaussianNumber(sqrt(T), 0.);
  }
}

void EvaluateTemperature (int sw, double pressure, double kin_energy, int *NT, double *t_eq, double *P) {
  if((*NT)==0) printf("Starting evaluating temperature at equilibrium...\n");

  switch (sw) {
    case 0:
      (*NT)++;
      (*t_eq) += (2.*kin_energy/DIM);
      (*P) += pressure;
      break;

    case 1:
      printf("\nSimulation completed!\n");
      *t_eq /= (*NT);
      *P /= (*NT);
      break;
  }
}

void Distances(double d[27], vec r1, vec r2) {
  int i, j, k, c = 0;
  double R2 = SquaredNorm(SubVec(r1, r2));
  for(i=-1; i<=1; i++) {
    for(j=-1; j<=1; j++) {
      for(k=-1; k<=1; k++) {
        d[c] = sqrt(R2 + doubleL*(i*(r1.x-r2.x) +  j*(r1.y-r2.y) + k*(r1.z-r2.z)) + L_2*(abs(i)+abs(j)+abs(k)));
        c++;
      }
    }
  }
}

void ModifyTemperature (double *T) {
  // if((*T)>2.) (*T) -= 0.001;
  // if((*T)<2. && (*T)>0.001) (*T) -= 0.001;
  if((*T)>0.01) (*T) -= 0.001;
  // *T = 4;
}

void EnergyAndPressureFromRDF (double g[], double *pressure, double *pot_energy, double t_eq) {
  int i;
  double r, p = 0., u = 0.;
  for(i=0; i<nBins_RDF; i++) {
    r = binsize_RDF*(i + 0.5);
    u += (g[i]*binsize_RDF*LJPotential(r)*r*r);
    p += (g[i]*binsize_RDF*r*r*r*DerivativeLJ(r));
  }
  u *= 2*M_PI*rho;
  p *= (-2./3*M_PI*rho*rho);
  p += rho*t_eq;
  *pressure = p;
  *pot_energy = p;
  printf("Potential Energy from g(r): %.10lf\n", u);
  // printf("Pressure from g(r): %lf\n", p);
}

void PrintVelocityField (mol particles[]) {
  int i;
  FILE *p = fopen("vfield.dat", "w");
  for (i=0; i<nParticles; i++) {
    fprintf(p, "%lf %lf %lf ", particles[i].r.x, particles[i].r.y, particles[i].r.z);
    fprintf(p, "%lf %lf %lf\n", particles[i].v.x, particles[i].v.y, particles[i].v.z);
  }
}

// void Density(mol particles[]) {
//   printf("%d\n", N_DEN);
//   double M[N_DEN][N_DEN][N_DEN];
//   int i, j, k;
//   for(i=0; i<N_DEN; i++) {
//     for(j=0; j<N_DEN; j++) {
//       for(k=0; k<N_DEN; k++) {
//         M[i][j][k] = 0;
//       }
//     }
//   }
//
//   vec p[nParticles];
//
//   int x, y, z;
//   FILE *d = fopen("density.dat", "w");
//
//   for(i=0; i<nParticles; i++) {
//     p[i] = SetEqualVec(particles[i].r);
//     p[i].x += 0.5*L;
//     p[i].y += 0.5*L;
//     p[i].z += 0.5*L;
//     x = (int)(p[i].x/binsize_DEN);
//     y = (int)(p[i].y/binsize_DEN);
//     z = (int)(p[i].z/binsize_DEN);
//     printf("%d %d %d\n", x, y, z);
//     M[x][y][z] += 1.;
//   }
//
//   double sum = 0.;
//
//   for(i=0; i<N_DEN; i++) {
//     for(j=0; j<N_DEN; j++) {
//       for(k=0; k<N_DEN; k++) {
//         fprintf(d, "%lf %lf %lf %lf\n", i*binsize_DEN-0.5*L, j*binsize_DEN-0.5*L, k*binsize_DEN-0.5*L, M[i][j][k]/(double)nParticles);
//         sum += M[i][j][k];
//       }
//       fprintf(d, "\n");
//     }
//     fprintf(d, "\n");
//   }
//
//   printf("%lf\n", sum);
//
// }

void Density(mol particles[]) {
  int i, j, x, y;
  vec p[nParticles];
  SetZeroArrayVec(p, nParticles);
  double M[N_DEN][N_DEN];

  FILE *d = fopen("d.dat", "w");

  for(i=0; i<N_DEN; i++) {
    for(j=0; j<N_DEN; j++) {
      M[i][j] = 0.;
    }
  }

  for(i=0; i<nParticles; i++) {
    p[i] = SetEqualVec(particles[i].r);
    p[i].x += 0.5*L;
    p[i].y += 0.5*L;
    p[i].z += 0.5*L;
    if(p[i].z>0.5*L && p[i].z<(0.5*L+0.5)) {
      x = (int)(p[i].x/binsize_DEN);
      y = (int)(p[i].y/binsize_DEN);
      // printf("%d %d\n", x, y);
      M[x][y] += 1.;
    }
  }

  for(i=0; i<N_DEN; i++) {
    for(j=0; j<N_DEN; j++) {
      M[i][j] /= (double)nParticles;
      fprintf(d, "%lf %lf %lf\n", (i)*binsize_DEN-0.5*L, (j)*binsize_DEN-0.5*L, M[i][j]);
    }
    fprintf(d, "\n");
  }

}
