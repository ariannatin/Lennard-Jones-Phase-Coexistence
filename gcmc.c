#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "parameters.h"
#include "MClib.h"
#include "vector_operations.h"
#include "sample.h"
// #include "physics.h"

int main() {
  // clock starts
  clock_t begin=clock();

  // initialization random seeds
  srand((long)time(NULL));
  srand48((long)time(NULL));

  // variables
  vec p[NMAX];                                          // position of every particle
  int i, k;                                             // step
  int nP = 0;                                           // number of particles
  int nPin;                                             // number of particles at the beginning
  int count_displ = 0, count_del = 0, count_ins = 0;    // counters of accepted trials
  int att_displ = 0, att_ex = 0;                        // counters of total trials

  // opening files
  FILE *pos = fopen("pos.dat", "w");
  FILE *num = fopen("num.dat", "w");
  FILE *en = fopen("en.dat", "w");
  FILE *init = fopen("init_pos.dat", "r");
  FILE *last_conf = fopen("last_conf.dat", "w");
  FILE *altro = fopen("400p.dat", "r");

  // initialization
  InitConfiguration(altro, p, &nP);
  //PrintAllPos(p, nP, 0);
  // PrintAllPosFile(pos, p, nP, 0);
  nPin = nP;

  // Monte Carlo steps
  for(i=0; i<nSteps; i++) {
    for(k=0; k<(nEx+nDispl); k++) {
      if(RandInt(1, (nEx+nDispl))<=nDispl){
        AttemptParticleDisplacement(p, nP, &count_displ);
        att_displ++;
      } else {
        AttemptParticleExchange(p, &nP, &count_del, &count_ins);
        // AttemptParticleInsertion(p, &nP, &count_ins);
        att_ex++;
      }
    }

    if(i%5 == 0) SampleNumberOfParticles(num, nP, i);

    if((i+1)%100 == 0) printf("%d steps done...   %d particles\n", i+1, nP);
  }

  // Monte Carlo - canonical ensemble
  // for(i=0; i<nSteps; i++) {
  //   for(k=0; k<(nEx+nDispl); k++) {
  //     AttemptParticleDisplacement(p, nP, &count_displ);
  //     att_displ++;
  //     PrintEnergyFile(en, i, p, nP);
  //   }
  //   if(i%20 == 0) {
  //     printf("%d steps done...\n", i);
  //     printf("\nenergia: %lf\n", TotEnergy(p, nP));
  //   }
  // }

  // frequency of acceptance
  printf("tot displ %d\n", att_displ);
  // printf("tot ex %d\n", att_ex);
  printf("d %lf\n", (double)count_displ/att_displ);
  // printf("ex %lf\n", (double)(count_ins+count_del)/att_ex);

  PrintAllPosFile(pos, p, nP, 0);

  fprintf(last_conf, "%d\n", nP);
  for(i=0; i<nP; i++) {
    fprintf(last_conf, "%.20lf %.20lf %.20lf\n", p[i].x, p[i].y, p[i].z);
  }

  // informations about the simulation
  clock_t end=clock();
  double ex_time=(double)(end-begin)/CLOCKS_PER_SEC;
  Info(nPin, nP, ex_time);

  return 0;
}
