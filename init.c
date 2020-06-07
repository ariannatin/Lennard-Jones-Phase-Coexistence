#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "parameters.h"
#include "MClib.h"
#include "vector_operations.h"
#include "sample.h"

int main() {
  clock_t begin=clock();

  srand((long)time(NULL));
  srand48((long)time(NULL));

  FILE *pos = fopen("init_pos.dat", "w");

  vec p[nMax];                       // position of every particle
  int nP = 120;
  int i;

  InitPosition(p, nP);               // max 250/500 particles
  // SCInitialization(p, &nP);       // 343/748 particles
  // BCCInitialization(p, &nP);      // 686/1495 particles
  // FCCInitialization(p, &nP);      // 1372/2990 particles

  fprintf(pos, "%d\n", nP);

  for(i=0; i<nP; i++) {
    fprintf(pos, "%.20lf %.20lf %.20lf\n", p[i].x, p[i].y, p[i].z);
  }
  fprintf(pos, "\n\n");



}
