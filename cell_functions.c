#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cell_functions.h"

int NumberOfCellFromPosition (vec r) {
  vec v = SumVec(DivVec(r, r_cut), SetVec((double) ((int)(W/2)), (double)((int)(W/2)), (double)((int)(W/2))));

 return 1 + (int)(v.x) + (int)(v.y)*W + (int)(v.z)*W*W;
}

int NumberOfCell (vec v) {
  return 1 + (int)(v.x) + (int)(v.y)*W + (int)(v.z)*W*W;
}

void CreateCellList (int cellnumber[], int head[], int after[], int before[], mol particles[]) {
  int i, icell;

  for(i=0; i<=(W*W*W); i++) head[i] = 0;

  for(i=1; i<=nParticles; i++) {
    icell = NumberOfCellFromPosition(particles[i-1].r);
    cellnumber[i] = icell;
    icell = cellnumber[i];
    after[i] = head[icell];
    before[head[icell]] = i;
    head[icell] = i;
  }
}

void UpdateCellParticle (int head[], int after[], int before[], int i, int old_icell, int new_icell) {
  int tmp;

  if(new_icell!=old_icell) {
    if(before[i]==0) head[old_icell] = after[i];
    // printf("%d \n", new_icell);
    tmp = head[new_icell];
    head[new_icell] = i;
    after[before[i]] = after[i];
    before[after[i]] = before[i];
    after[i] = tmp;
    before[i] = 0;
    before[after[i]] = i;
  }
}

void UpdateCellList (int head[], int after[], int before[], int cellnumber[], mol particles[]) {
  int i, new_icell, old_icell;

    for (i=1; i<=nParticles; i++) {
      new_icell = NumberOfCellFromPosition(particles[i-1].r);
      if(new_icell<0) printf("!!!!!!!!!!!!!\n");
      // PrintVec(particles[i-1].r);
      // printf("particella %d, newicell %d \n", i, new_icell);
      // PrintVec(particles[i-1].r);
      old_icell = cellnumber[i];
      // printf("%d %d %d\n", i, old_icell, new_icell);
      cellnumber[i] = new_icell;
      UpdateCellParticle(head, after, before, i, old_icell, new_icell);
  }
}
