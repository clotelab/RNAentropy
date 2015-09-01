#include "myConst.h"

double fTemp; // Formal temperature
double fkT;   // RT of formal temperature

int McGetZB(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZM, double **ZB);
int McGetZM1(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZB);
int McGetZM(int i, int j, char sequence[MAXSIZE],double **ZM1, double **ZM);
int McGetZ(int i, int j, char sequence[MAXSIZE],double **Z,double **ZB);
double** runMcCaskill(char sequence[MAXSIZE]);

int McGetQB(int i,int j, char sequence[MAXSIZE], double **QM1, double **QM, double **QB, double **ZM1, double **ZM, double **ZB);
int McGetQM1(int i, int j, char sequence[MAXSIZE], double **QM1, double **QB, double **ZM1, double **ZB);
int McGetQM(int i, int j, char sequence[MAXSIZE], double **QM1, double **QM, double **ZM1, double **ZM);
int McGetQ(int i, int j, char sequence[MAXSIZE], double **Q,double **QB, double **Z,double **ZB);

void mcCaskillQ(char sequence[MAXSIZE], double **Q, double **QB, double **QM1, double **QM, double **Z, double **ZB, double **ZM1, double **ZM);
void mcCaskill(char sequence[MAXSIZE], double **Z, double **ZB, double **ZM1, double **ZM);

int basePair(int i,int j,char rna[MAXSIZE]);
