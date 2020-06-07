#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "MClib.h"

void AttemptParticleDisplacement (vec p[], int nP, int *accepted) {
  int o;
  double e_old, e_new, a, b, c;
  vec prop = SetVec(0., 0., 0.);

  // random selection of a particle
  o = RandInt(0, nP-1);
  // PrintI(o, "o");
  // energy of the old configuration
  e_old = Energy(p, nP, p[o], o);

  // random displacement of particle
  a = RandReal(-0.5, 0.5);
  b = RandReal(-0.5, 0.5);
  c = RandReal(-0.5, 0.5);

  // PrintVecAndName(p[o], "prima");
  prop = SumVec(p[o], MultiplyVec(SetVec(a, b, c), delta));
  // PrintVecAndName(prop, "dopo 1");
  prop = PeriodicImage(prop);
  // PrintVecAndName(prop, "dopo 2");

  // energy of the new configuration
  e_new = Energy(p, nP, prop, o);

  // accept or reject?
  if(e_new<e_old) {
    // PrintD(e_old-e_new, "diff 1");
    p[o] = prop;
    (*accepted) ++;
  }

  else if(RandReal(0, 1) < exp(beta*(e_old-e_new))) {
    p[o] = prop;
    (*accepted) ++;
    // PrintD(exp(beta*(e_old-e_new)), "diff 2");
  }

}

void AttemptParticleDeletion (vec p[], int *nP, int *accepted) {
  if(*nP != 0) {
    double e_old, arg;
    int o;

    // random selection of a particle
    o = RandInt(0, *nP -1);

    // energy of the old configuration
    e_old = Energy(p, *nP, p[o], o);

    // acceptance rule
    arg = (*nP) * exp(beta*e_old) / (zz*V);
    // printf("\narg: %lf, %lf e_old, zz %lf\n", arg, e_old, zz);

    // PrintAllPos(p, *nP, 0);
    // Space();
    // if accepted, remove the particle
    if(RandReal(0, 1)<arg) {
      p[o] = p[*nP -1];
      p[*nP -1] = SetVec(0., 0., 0.);
      *nP = (*nP)-1;
      // printf("Accepted rimozione!\n");
      // PrintAllPos(p, *nP, 0);
      // Space();
      (*accepted) ++;
    }
  }
}

void AttemptParticleInsertion (vec p[], int *nP, int *accepted) {
  double a, b, c, e_new, arg;

  // new particle at a random position
  if(*nP<(NMAX-1)) {
    a = RandReal(-L/2., L/2.);
    b = RandReal(-L/2., L/2.);
    c = RandReal(-L/2., L/2.);
    vec prop = SetVec(a, b, c);

    // energy new configuration
    e_new = Energy(p, *nP, prop, (*nP+10));

    // acceptance rule
    arg = zz * V * exp(-beta*e_new) / ((*nP)+1);
    // PrintD(e_new, "e_new");
    // printf("\narg: %lf\n", arg);

    // if accepted, add the particle (in the last position)
    if (RandReal(0, 1)<arg) {
      // printf("\narg: %lf\n", arg);
      // PrintD(e_new, "e_new");
      // PrintVecAndName(prop, "prop");
      p[(*nP)] = prop;
      *nP = (*nP)+1;
      // printf("accettato aggiunta! \n");
      (*accepted) ++;
    }
  }
}

void AttemptParticleInsertionUS (vec p[], int *nP, int *accepted, int Nmin, int Nmax) {
  double a, b, c, e_new, arg;

  // new particle at a random position
  if(*nP<(NMAX-1) && *nP<Nmax) {
    a = RandReal(-L/2., L/2.);
    b = RandReal(-L/2., L/2.);
    c = RandReal(-L/2., L/2.);
    vec prop = SetVec(a, b, c);

    // energy new configuration
    e_new = Energy(p, *nP, prop, (*nP+10));

    // acceptance rule
    arg = zz * V * exp(-beta*e_new) / ((*nP)+1);
    // PrintD(e_new, "e_new");
    // printf("\narg: %lf\n", arg);

    // if accepted, add the particle (in the last position)
    if (RandReal(0, 1)<arg) {
      // printf("\narg: %lf\n", arg);
      // PrintD(e_new, "e_new");
      // PrintVecAndName(prop, "prop");
      p[(*nP)] = prop;
      *nP = (*nP)+1;
      // printf("accettato aggiunta! \n");
      (*accepted) ++;
    }
  }
}

void AttemptParticleDeletionUS (vec p[], int *nP, int *accepted, int Nmin, int Nmax) {
  if(*nP != 0 && *nP>Nmin) {
    double e_old, arg;
    int o;

    // random selection of a particle
    o = RandInt(0, *nP -1);

    // energy of the old configuration
    e_old = Energy(p, *nP, p[o], o);

    // acceptance rule
    arg = (*nP) * exp(beta*e_old) / (zz*V);
    // printf("\narg: %lf, %lf e_old, zz %lf\n", arg, e_old, zz);

    // PrintAllPos(p, *nP, 0);
    // Space();
    // if accepted, remove the particle
    if(RandReal(0, 1)<arg) {
      p[o] = p[*nP -1];
      p[*nP -1] = SetVec(0., 0., 0.);
      *nP = (*nP)-1;
      // printf("Accepted rimozione!\n");
      // PrintAllPos(p, *nP, 0);
      // Space();
      (*accepted) ++;
    }
  }
}

void AttemptParticleExchange (vec p[], int *nP, int *count_del, int *count_ins) {
  if(RandReal(0, 1)<0.5) {
    AttemptParticleDeletion (p, nP, count_del);
    // printf("removal\n");
  }
  else {
    AttemptParticleInsertion (p, nP, count_ins);
    // printf("add\n");
  }
}

void AttemptParticleExchangeUS (vec p[], int *nP, int *count_del, int *count_ins, int Nmin, int Nmax) {
  if(RandReal(0, 1)<0.5) {
    AttemptParticleDeletionUS (p, nP, count_del, Nmin, Nmax);
  }
  else {
    AttemptParticleInsertionUS (p, nP, count_ins, Nmin, Nmax);
  }
}


void FCCInitialization(vec p[], int *nP) {
  int i, j, k, l;
  int ampl = 1; //sigma
  int n = L;
  int number=0;
  int d = n/2;

  vec a[4];
  a[0] = SetVec(0., 0., 0.);
  a[1] = SetVec(ampl/2., ampl/2., 0.);
  a[2] = SetVec(0., ampl/2., ampl/2.);
  a[3] = SetVec(ampl/2., 0., ampl/2.);

  for(k=0; k<n; k++) {
    for(j=0; j<n; j++) {
      for(i=0; i<n; i++){
        vec v = MultiplyVec(SetVec(i, j, k), ampl);
        // PrintVec(v);
        for(l=0; l<4; l++){
          p[number] = SumVec(v, a[l]);
          p[number] = SubVec(p[number], SetConstantVector(L/2.));
          p[number] = PeriodicImage(p[number]);
          number++;
       }
      }
    }
  }
  //
  // i=7;
  // for(j=0; j<n; j++) {
  //   for(k=0; k<n; k++) {
  //     vec v = MultiplyVec(SetVec(i, j, k), ampl);
  //     p[number] = SumVec(v, a[0]);
  //     number++;
  //
  //     p[number] = SumVec(v, a[2]);
  //     number++;
  //   }
  // }
  //
  // j=7;
  // for(i=0; i<n; i++) {
  //   for(k=0; k<n; k++) {
  //     vec v = MultiplyVec(SetVec(i, j, k), ampl);
  //     p[number] = SumVec(v, a[0]);
  //     number++;
  //
  //     p[number] = SumVec(v, a[3]);
  //     number++;
  //   }
  // }
  //
  // k=7;
  // for(j=0; j<n; j++) {
  //   for(i=0; i<n; i++) {
  //     vec v = MultiplyVec(SetVec(i, j, k), ampl);
  //     p[number] = SumVec(v, a[0]);
  //     number++;
  //
  //     p[number] = SumVec(v, a[1]);
  //     number++;
  //   }
  // }
  //
  // k=7;
  // i=7;
  // for(j=0; j<n; j++){
  //   vec v = MultiplyVec(SetVec(i, j, k), ampl);
  //   p[number] = v;
  //   number++;
  // }
  //
  // k=7;
  // j=7;
  // for(i=0; i<n; i++) {
  //   vec v = MultiplyVec(SetVec(i, j, k), ampl);
  //
  //   p[number] = v;
  //
  //   number++;
  // }
  //
  // i=7;
  // j=7;
  // for(k=0; k<n; k++) {
  //   vec v = MultiplyVec(SetVec(i, j, k), ampl);
  //
  //   p[number] = v;
  //
  //   number++;
  // }
  //
  // i=7;
  // j=7;
  // k=7;
  // p[number] = MultiplyVec(SetVec(i, j, k), ampl);
  // number++;

  *nP = number;

}

void InitRandomPosition (vec p[], int nP) {
  int k = 0, n = 0, ll = 0;
  vec tmp;
  double l;

  printf("\n\nParticles initialized:\n0\n");

  while (n != nP) {
    ll = 0;
    tmp = SetVec(RandReal(-L/2., L/2.), RandReal(-L/2., L/2.), RandReal(-L/2., L/2.));
    for(k=0; k<n; k++) {
      l = SquaredNorm(SeparationPBC(tmp, p[k]));
      if (l<1.) {
        ll = 1;
        break;
      }
    }
    if(ll!=1) {
      p[n] = SetEqualVec(tmp);
      n++;
      if(n%10 == 0) printf("%d \n", n);
    }
  }
  printf("===========================\n\n\nNow the simulation begins:\n\n\n");
  //
  // for(k=0; k<nP; k++) {
  //   p[k] = DivVec(p[k], L);
  // }
}

void SCInitialization(vec p[], int *nP) {
  int i, j, k;
  int ampl = 1; //sigma
  int n = L;
  int number=0;
  int d = n/2;
  vec v;

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      for(k=0; k<n; k++) {
        p[number] = MultiplyVec(SetVec(i, j, k), ampl);
        p[number] = SubVec(p[number], SetConstantVector(L/2.));
        p[number] = PeriodicImage(p[number]);
        number++;
      }
    }

    *nP = number;
  }
}

void BCCInitialization(vec p[], int *nP) {
  int i, j, k, l;
  int ampl = 1; //sigma
  int n = L;
  int number=0;
  int d = n/2;

  vec a[4];
  a[0] = SetVec(0., 0., 0.);
  a[1] = SetVec(ampl/2., ampl/2., ampl/2.);

  for(k=0; k<n; k++) {
    for(j=0; j<n; j++) {
      for(i=0; i<n; i++){
        vec v = MultiplyVec(SetVec(i, j, k), ampl);
        for(l=0; l<2; l++){
          p[number] = SumVec(v, a[l]);
          p[number] = SubVec(p[number], SetConstantVector(L/2.));
          p[number] = PeriodicImage(p[number]);
          number++;
        }
      }
    }
  }

  *nP = number;

}

void InitConfiguration (FILE*point, vec p[], int *nP) {
  int i;
  double a, b, c;

  if(INIT_FROM_FILE == 1) {
    // fscanf(point, "%d", nP);
    if(fscanf(point, "%d", nP)){
      printf("%d\n", *nP);
      for(i=0; i<(*nP); i++) {
        if(fscanf(point, "%lf %lf %lf", &a, &b, &c)){
          p[i] = SetVec(a, b, c);
        }
      }
    }
  }
  else {
      InitRandomPosition(p, *nP);
  }
}


void HistogramOfN (double H[], double R[], double HR[], double P[], FILE *prob){
  int i, n;
  double sum = 0;
  double Pmax = -1000;
  FILE *problog = fopen("log.dat", "w");
  SetZeroRealArray(P, N_SAMPLE_MAX+1);

  for(i=0; i<N_SAMPLE_MAX+1; i++) {
    if(H[(int)(i/W) * W] == 0) H[(int)(i/W) * W] = 0.00000001;
    HR[i] = H[i]/(double)H[(int)(i/W) * W];
  }

  PrintRealArray("HR", HR, N_SAMPLE_MAX);

  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    for(i=0; i<(int)(n/W); i++) {
      P[n] += log(R[i]);
    }
    P[n] += log(HR[n]);
    if(P[n]>Pmax) Pmax = P[n];
  }

  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    P[n] -= Pmax;
    P[n] = exp(P[n]);
    sum += P[n];
  }

  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    P[n] /= sum;
    fprintf(prob, "%d %.5e\n", n, P[n]);
    fprintf(problog, "%d %.5e\n", n, log(P[n]));
  }

}

void HistogramOfN_rec (double H[], double F[], double P[], FILE *prob) {
  int i, n, f = 0, index;
  double sum = 0;
  double Pmax = -1000;
  double HR[N_SAMPLE_MAX+1];
  double R[N_SAMPLE_MAX+1];
  SetZeroRealArray(P, N_SAMPLE_MAX+1);
  SetZeroRealArray(HR, N_SAMPLE_MAX+1);
  SetZeroRealArray(R, N_WINDOWS);

  for(i=0; i<N_WINDOWS; i++) {
    if(F[f] == 0) {
      F[f] = 0.00000001;
      printf("Attenzione!\n");
    }
    R[i] = F[f+1]/F[f];
    // printf("%lf %lf\n", F[f+1], F[f]);
    f+=2;
  }

  // PrintRealArray("R", R, N_WINDOWS);

  for(i=0; i<N_SAMPLE_MAX+1; i++) {
    index = (int)(i/W) * W;
    // printf("index %d\n", index);
    if(H[index] == 0) {
      H[index] = 0.00000001;
      // printf("Attenzione!\n");
    }
    // printf("%lf %lf\n", H[i], H[index]);
    HR[i] = H[i]/(double)H[index];
  }

  // PrintRealArray("HR", HR, N_SAMPLE_MAX+1);

  P[0] = 1;
  for(n=1; n<N_SAMPLE_MAX+1; n++) {
    for(i=0; i<(int)(n/W); i++) {
      P[n] += log(R[i]);
      // printf("n %d i %d\n", n, i);
    }
    P[n] += log(HR[n]);
    if(P[n]>Pmax) Pmax = P[n];
  }

  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    P[n] -= Pmax;
    P[n] = exp(P[n]);
    sum += P[n];
  }

  printf("somma: %lf\n", sum);

  double SUM = 0;

  for(n=0; n<N_SAMPLE_MAX+1; n++) {
    P[n] /= sum;
    SUM += P[n];
    fprintf(prob, "%d %.5e\n", n, P[n]);
  }

  printf(" SUM %lf\n", SUM);

  // PrintRealArray("P", P, N_SAMPLE_MAX+1);
}


void SaveLastSimulation (double H[], double F[], int run) {
  int m;

  char histo[50];
  sprintf(histo, "dati/histo/h%d.dat", run);
  FILE *h = fopen(histo, "w");
  printf("genero %s\n\n", histo);
  for(m=0; m<(N_SAMPLE_MAX+1); m++) {
    fprintf(h, "%e\n", H[m]);
  }

  char fr[50];
  sprintf(fr, "dati/fract/f%d.dat", run);
  FILE *fra = fopen(fr, "w");
  printf("genero %s\n\n", fr);
  for(m=0; m<(2*N_WINDOWS); m++) {
    fprintf(fra, "%e\n", F[m]);
  }
}

void SaveLastConfiguration (vec p[], int nP, int Nmin, int Nmax, int run) {
  char config[50];
  sprintf(config, "dati/conf/win%d_%d-%d.dat", run, Nmin, Nmax);
  FILE *conf = fopen(config, "w");
  printf("genero %s\n\n", config);
  // PrintAllPos(p, nP, 0);
  fprintf(conf, "%d\n", nP);
  PrintAllPosFile(conf, p, nP, 0);
}

void UploadLastConfiguration (vec p[], int nP, int Nmin, int Nmax, int run) {
  double a, b, c;
  int l;

  if(run == 0) return;

  char config[50];
  sprintf(config, "dati/conf/win%d_%d-%d.dat", run-1, Nmin, Nmax);
  FILE *conf = fopen(config, "r");
  printf("leggo da %s\n", config);
  if(fscanf(conf, "%d", &nP));
  // printf("nP: %d\n", nP);
  for(l=0; l<nP; l++) {
    if(fscanf(conf, "%lf %lf %lf", &a, &b, &c)){
      p[l] = SetVec(a, b, c);
    }
  }
}

void UploadLastSimulation (double H[], double F[], int run) {
  int m;

  if(run == 0) return;

  char histo[50];
  sprintf(histo, "dati/histo/h%d.dat", run-1);
  FILE *h = fopen(histo, "r");
  printf("leggo da %s\n", histo);
  for(m=0; m<(N_SAMPLE_MAX+1); m++) {
    if(fscanf(h, "%lf", &H[m]));
  }

  char fr[50];
  sprintf(fr, "dati/fract/f%d.dat", run-1);
  FILE *fra = fopen(fr, "r");
  printf("leggo da %s\n", fr);
  for(m=0; m<(2*N_WINDOWS); m++) {
    if(fscanf(fra, "%lf", &F[m]));
  }

  Space();
}
