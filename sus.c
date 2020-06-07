#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "parameters.h"
#include "MClib.h"
#include "vector_operations.h"
#include "sample.h"

int main(){
  // clock starts
  clock_t begin=clock();

  // initialization random seeds
  srand((long)time(NULL));
  srand48((long)time(NULL));

  // variables
  vec p[NMAX];                                          // position of every particle
  int i, k, j;                                          // step
  int n;                                                // particle index
  int nP = 0;                                           // number of particles
  int nPin;                                             // number of particles at the beginning
  int count_displ = 0, count_del = 0, count_ins = 0;    // counters of accepted trials
  int att_displ = 0, att_ex = 0;                        // counters of total trials
  double Hr = 0., Hl = 0.;                              // number of particles counters
  int Nmin = 10, Nmax = 20;                             // minimum and maximum number of particles for a window
  double sum = 0.;                                      // normalization of P(N)

  double R[N_WINDOWS], HR[N_SAMPLE_MAX], P[N_SAMPLE_MAX], H[N_SAMPLE_MAX];
  SetZeroRealArray(R, N_WINDOWS);
  SetZeroRealArray(HR, N_SAMPLE_MAX);
  SetZeroRealArray(P, N_SAMPLE_MAX);
  SetZeroRealArray(H, N_SAMPLE_MAX);

  // files
  FILE *pos = fopen("posUS15.dat", "w");
  FILE *num = fopen("numUS15.dat", "w");
  FILE *altro = fopen("400p15.dat", "r");
  FILE *prob = fopen("pdin15.dat", "w");

  FILE *p1 = fopen("p1.dat", "w");
  FILE *p2 = fopen("p2.dat", "w");
  FILE *p3 = fopen("p3.dat", "w");


  // initialization and estimated execution time
  double est_time = (12.0*nSteps/300.*W/2.*pow(N_SAMPLE_MAX/450., 3.5));
  printf("\n\n\n==========================\nEstimated execution time: \n%.3lf s = %.2lf min = %.2lf h\n", est_time*60, est_time, est_time/60.);

  InitConfiguration(altro, p, &nP);
  nPin = nP;

  // successive umbrella sampling
  for(j=0; j<N_WINDOWS; j+=1) {

    Hr = 0.;
    Hl = 0.;
    Nmin = j*W;
    Nmax = (j+1)*W;
    H[Nmin] = 0.;

    printf("number of particles: %d-%d\n", Nmin, Nmax);

    // Monte Carlo steps
    for(i=0; i<nSteps; i++) {
      for(k=0; k<(nEx+nDispl); k++) {
        if(RandInt(1, (nEx+nDispl))<=nDispl){
          AttemptParticleDisplacement(p, nP, &count_displ);
          att_displ++;
        } else {
          AttemptParticleExchangeUS(p, &nP, &count_del, &count_ins, Nmin, Nmax);
          att_ex++;
        }
      }

      if(nP == Nmax) Hr++;
      if(nP == Nmin) Hl++;

      H[nP]++;

      // fprintf(num, "%d\n", nP);
    }

    if(Hl == 0) Hl = 0.00000001;
    if(Hr == 0) Hr = 0.00000001;
    R[j] = Hr/(double)Hl;
  }

  PrintRealArray("H", H, N_SAMPLE_MAX);
  PrintRealArray("R", R, N_WINDOWS);

  HistogramOfN(H, R, HR, P, p1);

  // informations about the simulation
  clock_t end=clock();
  double ex_time=(double)(end-begin)/CLOCKS_PER_SEC;
  Info(nPin, nP, ex_time, P);

  return 0;
}
