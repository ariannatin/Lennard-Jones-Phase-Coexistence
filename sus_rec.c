#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "parameters.h"
#include "MClib.h"
#include "vector_operations.h"
#include "sample.h"

int main(int argc, char const *argv[]){
  // clock starts
  clock_t begin = clock();

  //  must decide, when executing, which run of the simulation
  if(argc == 2) {
    printf("\n\n\nRun %s of simulation\n", argv[1]);
  } else {
    printf("Insert the number of the run\n");
    return 0;
  }

  int run = atoi(argv[1]);

  // initialization random seeds
  srand((long)time(NULL));
  srand48((long)time(NULL));

  // =========================================================================== //
  // declaration of variables
  vec p[NMAX];                                          // position of every particle
  int i, k, j;                                          // step
  int n;                                                // particle index
  int nP = 0;                                           // number of particles
  int nPin;                                             // number of particles at the beginning
  int count_displ = 0, count_del = 0, count_ins = 0;    // counters of accepted trials
  int att_displ = 0, att_ex = 0;                        // counters of total trials
  int f = 0, l, m;                                      // auxiliary indices
  double a, b, c;                                       // auxiliary variables
  double Hr = 0., Hl = 0.;                              // number of particles counters
  int Nmin, Nmax;                                       // minimum and maximum number of particles for a window
  double sum = 0.;                                      // normalization of P(N)

  double P[N_SAMPLE_MAX+1], H[N_SAMPLE_MAX+1], F[2*N_WINDOWS];
  SetZeroRealArray(F, 2*N_WINDOWS);
  SetZeroRealArray(P, N_SAMPLE_MAX+1);
  SetZeroRealArray(H, N_SAMPLE_MAX+1);

  // opening files
  char pdn[30];
  sprintf(pdn, "dati/prob/pdn%d.dat", run);
  FILE *prob = fopen(pdn, "w");
  char free[30];
  sprintf(free, "dati/en/en%d.dat", run);
  FILE *en = fopen(free, "w");

  // initialization and estimated execution time
  double est_time = (12.0*nSteps/300.*3*W/2.*pow(N_SAMPLE_MAX/450., 3.5));
  printf("\n\n============================\nEstimated execution time: \n%.3lf s = %.2lf min = %.2lf h\n", est_time*60, est_time, est_time/60.);
  printf("============================\n\n\n");

  // the simulation starts with no particles
  nPin = 0;

  // =========================================================================== //

  // upload the selected run reading from the selected simulation
  // here the arrays that will be needed to estimate the density probability are loaded

  UploadLastSimulation(H, F, run);
  // PrintRealArray("H", H, N_SAMPLE_MAX+1);
  // PrintRealArray("F", F, 2*N_WINDOWS);

  double tmpH = H[0];

  // successive umbrella sampling:
  // a simulation is performed for every window with a fixed range of particles
  for(j=0; j<N_WINDOWS; j+=1) {

    Hr = 0.;
    Hl = 0.;
    Nmin = j*W;
    Nmax = (j+1)*W;

    H[Nmin] = tmpH;
    tmpH = H[Nmax];

    printf("number of particles: %d-%d\n", Nmin, Nmax);

    // initialization
    // here the last configuration of the selected run is loaded for every window
    UploadLastConfiguration(p, nP, Nmin, Nmax, run);

    // Monte Carlo steps for Grand Canonical Ensemble in [Nmin, Nmax] particles
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

    }

    F[f] += Hl;
    f++;

    F[f] += Hr;
    f++;

    // save configurations for every window
    SaveLastConfiguration(p, nP, Nmin, Nmax, run);

    // PrintRealArray("H", H, N_SAMPLE_MAX+1);
    // PrintRealArray("F", F, 2*N_WINDOWS);
  }

  // =========================================================================== //

  // save the arrays for the estimation of the probability density
  SaveLastSimulation(H, F, run);

  // compute histogram of the number of particles
  HistogramOfN_rec(H, F, P, prob);

  // compute the free energy
  FreeEnergy(P, en);

  // histogram reweighting in zeta
  FILE *newp = fopen("altrap.dat", "w");
  double newP[N_SAMPLE_MAX+1];
  double new_zeta = 0.047;
  HistogramReweightingZ(new_zeta, P, newP);
  // PrintRealArray("newp", newP, N_SAMPLE_MAX+1);
  for(i=0; i<N_SAMPLE_MAX+1; i++) {
    fprintf(newp, "%lf\n", newP[i]);
  }

  // =========================================================================== //

  // print informations about the simulation
  clock_t end = clock();
  double ex_time = (double)(end-begin)/CLOCKS_PER_SEC;
  Info(nPin, nP, ex_time, P);

  return 0;
}
