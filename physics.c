#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "physics.h"

void Info (mol particles[]) {
  printf("\n\n*******************************************************\n");
  printf("  This is a simulation of a Lennard-Jones fluid\n  with the following parameters:");
  printf("\n=======================================================\n");
  printf("       Particles: %d          \n", nParticles);
  printf("       Density: %lf                 \n", rho);
  printf("       Cube of %lf x %lf x %lf                 \n", L, L, L);
  printf("       Number of cells per side: %d            \n", W);
  printf("       Initial temperature %lf                \n", TEMPERATURE);
  printf("       Absolute temperature %lf                \n", TT);
  printf("       Integration stops at %lf               \n", TMAX);
  printf("       Number of steps: %lf                  \n", (double)TMAX/dt);
  printf("=======================================================\n");
}

double LJPotential(double r) {
  double r2 = r*r;
  if(r2<r2_cut) {
    double r2i = 1./r2, r6i = Cube(r2i);
    return 4*r6i*(r6i - 1.) + 4*e_cut;
  } else return 0.;
}

double DerivativeLJ(double r) {
  double r2 = r*r;
  if(r2<r2_cut) {
    double r2i = 1./r2, r6i = Cube(r2i), r12i = r6i*r6i;
    return 24./r*(-2.*r12i + r6i);
  } else return 0.;
}

void Neighbours (int neigh[1+nCells][14], vec cell) {
  vec neigh_temp;
  double i, j, k;
  int n = 0, icell = NumberOfCell(cell);

  for(i=-1; i<=0; i+=1) {
    for(j=-1; j<=1; j+=1){
      for(k=-1; k<=1; k+=1){

        if  (!(j==-1 && i==0) && !(i==0 && j==0 && k==1)) {
          neigh_temp = SetEqualVec(cell);

          neigh_temp.x += i;
          neigh_temp.y += j;
          neigh_temp.z += k;

          if(neigh_temp.x >= W) neigh_temp.x -= W;
          if(neigh_temp.x < 0) neigh_temp.x += W;
          if(neigh_temp.y >= W) neigh_temp.y -= W;
          if(neigh_temp.y < 0) neigh_temp.y += W;
          if(neigh_temp.z >= W) neigh_temp.z -= W;
          if(neigh_temp.z < 0) neigh_temp.z += W;

          neigh[icell][n] = NumberOfCell(neigh_temp);

          // printf("icell %d neigh %d\n", icell, neigh[icell][n]);

          n++;
        }
      }
    }
  }
}

void BuildNeighbours (int neigh[1+nCells][14]) {
   vec cell;
    for(cell.x=0; cell.x<W; cell.x+=1) {
      for(cell.y=0; cell.y<W; cell.y+=1) {
        for(cell.z=0; cell.z<W; cell.z+=1) {
         Neighbours(neigh, cell);
        }
      }
    }
}

void Setup (mol particles[], double *kin_energy, double *pot_energy) {
  printf("Initializing positions...\n");
  InitPosition(particles);
  printf("Done!\n\n");
  InitVel(particles);
  (*kin_energy) = KinEnergy(particles); //mettere dentro setup
  (*pot_energy) = PotEnergy(particles);
}

void InitPosition (mol particles[]) {
  int k = 0, n = 0, ll = 0;
  vec tmp;
  double l;

  printf("Particles initialized: \n");

  while (n != nParticles) {
    ll = 0;
    tmp = SetVec(RandReal(-L/2., L/2.), RandReal(-L/2., L/2.), RandReal(-L/2., L/2.));
    // tmp = SetVec(0., RandReal(-L/2., L/2.), RandReal(-L/2., L/2.));
    for(k=0; k<n; k++) {
      l = SquaredNorm(SeparationPBC(tmp, particles[k].r));
      if (l<1.) {
        ll = 1;
        break;
      }
    }
    if(ll!=1) {
      particles[n].r = SetEqualVec(tmp);
      n++;
      if(n%10 == 0) printf("%d \n", n);
    }
  }
}

void InitVel (mol particles[]) {
  int i;
  vec meanV = SetVec(0, 0, 0), M;

  double meanV2=0, scalingFactor;

  for(i=0; i<nParticles; i++)  {
    particles[i].v = SetVec(RandReal(-1., 1.), RandReal(-1., 1.), RandReal(-1., 1.));
    // particles[i].v = SetVec(0., RandReal(-1., 1.), RandReal(-1., 1.));
    meanV = SumVec(meanV, (particles[i].v));
    meanV2 += SquaredNorm(particles[i].v);
  }

  meanV = DivVec (meanV, nParticles);
  meanV2 = meanV2/nParticles;
  scalingFactor = sqrt(3*TEMPERATURE/meanV2);

  M = SetVec (meanV.x, meanV.y, meanV.z);

  meanV = SetVec(0, 0, 0);
  meanV2 = 0;

  for(i=0; i<nParticles; i++)  {
    particles[i].v = SubVec((particles[i].v), M);
    particles[i].v = MultiplyVec((particles[i].v), scalingFactor);
    meanV = SumVec(meanV, (particles[i].v));
    meanV2 += SquaredNorm(particles[i].v);
  }

  meanV = MultiplyVec (meanV, 1./nParticles);
  meanV2 = meanV2/(2.*nParticles);
}

void FirstLeapfrogStep (mol particles[], int cellnumber[], int neigh[1+nCells][14], int head[1+nCells], int after[nParticles+1], int before[nParticles+1], double *pressure) {
  int i, k, l, j, jcell, ncell;
  vec forces[nParticles], ff;
  double pot_energy_particle = 0., p = 0., kin_energy = 0.;

  SetZeroArrayVec(forces, nParticles);

  printf("Starting the simulation...\n");
  for(i=0; i<nParticles; i++) {
    l = cellnumber[i+1];
    for(ncell=1; ncell<=14; ncell++) {
      if(ncell!=11) {
        jcell = neigh[l][ncell-1];
        j = head[jcell];
        while(j!=0) {
          if(i!=(j-1)) {
            ff = Force(particles[i].r, particles[j-1].r, &pot_energy_particle, &p);
            *pressure += p;
            forces[i] = SumVec(forces[i], ff);
            forces[j-1] = SubVec(forces[j-1], ff);
          }
          j = after[j];
        }
      }
    }
  }

  for(k=1; k<=nCells; k++) {
    i = head[k];
    while(i!=0) {
      j = after[i];
      while (j!=0) {
        ff = Force(particles[i-1].r, particles[j-1].r, &pot_energy_particle, &p);
        *pressure += p;
        forces[i-1] = SumVec(forces[i-1], ff);
        forces[j-1] = SubVec(forces[j-1], ff);
        j = after[j];
      }
      i = after[i];
    }
  }

  for (i=0; i<nParticles; i++) {
    kin_energy += SquaredNorm(particles[i].v);
    particles[i].v = SubVec(particles[i].v, MultiplyVec(forces[i], 0.5*dt));
  }

  kin_energy /= (2.*nParticles);
  *pressure /= DIM;
  *pressure += (nParticles*2.*(kin_energy)/DIM);
  *pressure /= Cube(L);
}

vec Force (vec r_i, vec r_j, double *pot_energy_couple, double *p) {

  vec f;;
  double r2, r2i, r6i;
  vec r_ij = SeparationPBC(r_i, r_j);
  r2 = SquaredNorm(r_ij);

  if(r2<r2_cut) {
    r2i = 1./r2;
    r6i = Cube(r2i);
    f = MultiplyVec(r_ij, (r6i * (r6i - 0.5) * r2i));
    *pot_energy_couple = r6i * (r6i - 1.) + e_cut;
    *p = ScalarProduct(r_ij, f);
  } else {
    f = SetVec(0., 0., 0.);
    *pot_energy_couple = 0.;
    *p = 0.;
  }

  if(r2<0.64) printf("Two particles got really close! Distance: %lf\n", sqrt(r2));
  return f;
}

void MoveParticle (mol particles[], vec displ[], vec f, int i, double T, double t, double EQ_TIME) {

  particles[i].v = SumVec(particles[i].v, MultiplyVec(f, 48. * dt));
  particles[i].r = SumVec(particles[i].r, MultiplyVec(particles[i].v, dt));

  particles[i].r = PeriodicImage(particles[i].r);
  if(t>EQ_TIME) displ[i] = SumVec(displ[i], MultiplyVec(particles[i].v, dt));

  #ifdef CONSTANT_TEMP
  AndersenThermostat(particles, i, T);
  #endif
}

void TakeStep (int cellnumber[], mol particles[], int head[], int after[], int neigh[1+nCells][14], double *pot_energy, double *kin_energy, vec* q, double *pressure, double *msd, vec displ[], double T, double t, double EQ_TIME) {
  int i = 0, j = 0, jcell = 0, ncell = 0, k = 0, l = 0;
  vec forces[nParticles], ff = SetVec(0., 0., 0.);
  double pot_energy_couple = 0., p = 0.;
  *kin_energy = 0.;
  *pot_energy = 0.;
  *msd = 0.;
  *pressure = 0.;

  SetZeroArrayVec(forces, nParticles);

  for(i=0; i<nParticles; i++) {
    l = cellnumber[i+1];
    for(ncell=1; ncell<=14; ncell++) {
      if(ncell!=11) {
        jcell = neigh[l][ncell-1];
        j = head[jcell];
        while(j!=0) {
          if(i!=(j-1)) {
            // printf("1\n");
            ff = Force(particles[i].r, particles[j-1].r, &pot_energy_couple, &p);
            *pot_energy += pot_energy_couple;
            *pressure += p;
            forces[i] = SumVec(forces[i], ff);
            forces[j-1] = SubVec(forces[j-1], ff);
            // printf("2\n");
          }
          j = after[j];
        }
      }
    }
  }

  for(k=1; k<=nCells; k++) {
    i = head[k];
    while(i!=0) {
        j = after[i];
        while (j!=0) {
        ff = Force(particles[i-1].r, particles[j-1].r, &pot_energy_couple, &p);
        *pot_energy += pot_energy_couple;
        *pressure += p;
        forces[i-1] = SumVec(forces[i-1], ff);
        forces[j-1] = SubVec(forces[j-1], ff);
        j = after[j];
      }
      i = after[i];
    }
  }

  vec sumv = {0, 0, 0};
  double m = 0.;
  for(i=0; i<nParticles; i++) {
    sumv = SumVec(sumv, particles[i].v);
    *kin_energy += SquaredNorm(particles[i].v);
    MoveParticle(particles, displ, forces[i], i, T, t, EQ_TIME);
    if(t>EQ_TIME)(*msd) += SquaredNorm(displ[i]);
  }

  *q = DivVec(sumv, nParticles);
  *pot_energy *= 4./(double)nParticles;
  *kin_energy /= (2.*nParticles);
  *pressure /= DIM;
  *pressure += (nParticles*2.*(*kin_energy)/DIM);
  *pressure /= Cube(L);
  if(t>EQ_TIME) *msd /= nParticles;
}

void AlternativeInitPosition (mol particles[]) {
  double i, j, k;
  int n = 0;
  int t = 0;

  // for(i=(-0.5*L); i<(0.5*L); i+=1.) {
  //   for(j=(-0.5*L); j<(0.5*L); j+=1.) {
  //     for(k=(-0.5*L); k<(0.5*L); k+=1.) {
  //       t++;
  //       if(RandReal(0., 1.)<rho) {
  //         particles[n].r = SetVec(i, j, k);
  //         n++;
  //       }
  //     }
  //   }
  // }

  // t = 0;
  //
  // for(i=0; i<(0.5*L); i+=1.) {
  //   for(j=0; j<(0.5*L); j+=1.) {
  //     for(k=0; k<(0.5*L); k+=1.){
  //       t++;
  //       printf("%lf\n", k);
  //     }
  //   }
  // }
  //
  // for(i=-1.; i>(-0.5*L); i-=1.) {
  //   for(j=-1.; j>(-0.5*L); j-=1.) {
  //     for(k=-1.; k>(-0.5*L); k-=1.){
  //       t++;
  //     }
  //   }
  // }
  t = 0;
  for(i=(-0.5*L); i<0.5*L; i+=1.) {
    printf("%lf\n", i);
    for(j=0; j<(0.5*L); j+=1.) {
      for(k=0; k<0.5*L; k+=1.) {
        t++;
      }
      for(k=-1.; k>-0.5*L; k-=1.){
        t++;
      }
    }
    for(j=-1.; j>-0.5*L; j-=1.) {
      for(k=0; k<0.5*L; k+=1.) {
        t++;
      }
      for(k=-1.; k>-0.5*L; k-=1.){
        t++;
      }
    }
  }

  printf("%d %lf\n", t, Cube(L));
  printf("%d %lf\n", n, n/Cube(L));
}
