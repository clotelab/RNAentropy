/* Program to calculate and store irred HP's */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "convert_Vienna.h"
#include "HP.h"
//#include "pair_mat.h"

#define EPSILON 0.00001
#define PRIVATE static
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

PRIVATE double **HP_energy_part;
PRIVATE double part, energy;
PRIVATE int is_min(int i, int j, short *S0, char *sequence, double kT);

void HP_init(int n){
  int i;
  HP_energy_part=(double**) malloc((n+1)*sizeof(double*));
  for (i=1;i<=n;i++)
    HP_energy_part[i]=(double *) calloc(n+1, sizeof(double));
}

void find_irred_HP(short *S0, char *sequence, double kT, int n){
  int i,j;
  //HP_energy_part=(double**) malloc((n+1)*sizeof(double*));
  //for (i=1;i<=n;i++)
  //  HP_energy_part[i]=(double *) calloc(n+1, sizeof(double));
  for (i=1;i<=n;i++)
    for (j=i+4;j<=n;j++)
      if (pair[S0[i]][S0[j]]){
	if (is_min(i,j, S0, sequence, kT)){
	  HP_energy_part[j][i]=part;
	}
	HP_energy_part[i][j]=energy;
      }
  return;
}
double HP_E(int i, int j){
  return HP_energy_part[i][j];
}
double HP_part(int i, int j){
  return HP_energy_part[j][i];
}
double HP_is_min(int i,int j){
  return HP_energy_part[j][i];
}
PRIVATE int is_min(int i, int j, short *S0, char *sequence, double kT){
  int ip, jp;
  
  energy=HP_Energy(i,j,S0, sequence);
  for (ip=i+1;ip<=MIN(i+30,j-1);ip++)
    for (jp=j-1;jp>=MAX(j-(30-(ip-i)),ip+4);jp--)
      if (pair[S0[ip]][S0[jp]]){
	if (IL_Energy(i,j,ip,jp, S0)+HP_Energy(ip,jp, S0,sequence)
	    +EPSILON<energy){
	  return 0;
	}
      }
  part=exp(-energy/kT);
  //printf("part = %f\n",part);
  //printf("ene = %f\n",energy);
  return 1;
}

  
