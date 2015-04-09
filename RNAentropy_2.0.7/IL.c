/* module that stores IL data */
#include <stdio.h>
#include <stdlib.h>
#include "IL.h"

#define PRIVATE static
#define SIZE_INIT 6

PRIVATE int n, i_last=0, j_last=0, index=0;
PRIVATE int **IL_table;
PRIVATE int pres_size=SIZE_INIT, size_inc=2;
PRIVATE void chk_alloc();

void IL_initialize(int n_in){
  int i;

  n=n_in;
  IL_table=(int**) malloc((n+1)*sizeof(int*));
  for (i=1;i<=n;i++)
    IL_table[i]=(int *) calloc(n+1, sizeof(int));
  IL = (struct IL_info *) malloc (pres_size*n*n*sizeof(struct IL_info));
  return;
}

void IL_reset(){
  int i,j;
  
  for (i=1;i<=n;i++)
    for (j=0;j<=n;j++)
      IL_table[i][j]=0;
}

void IL_add(double part, double E0, int i, int j, int ip, int jp, double count){
  if (i==i_last && j==j_last){
    index++; IL_table[j][i]++;
  } else {
    index++;IL_table[i][j]=index;IL_table[j][i]=1;
    i_last=i;j_last=j;
  }
  chk_alloc();
  IL[index].ip=ip;IL[index].jp=jp;IL[index].part=part;
  IL[index].energy=E0;IL[index].count=count;
  return;
}  
int IL_start(int i, int j){
  return IL_table[i][j];
}
int IL_end(int i, int j){
  return IL_table[i][j]+IL_table[j][i]-1;
}
int IL_index(){
  return index;
}

PRIVATE void chk_alloc(){
  if (index>=pres_size*n*n){
    pres_size+=size_inc;
    if (realloc(IL, pres_size*n*n*sizeof(struct IL_info))==NULL){
      printf("failed memory allocation.\n");
      exit(EXIT_FAILURE);
    }
  }
  return;
}
