/*
   prototypes for energy_par.c
*/

#ifndef __VIENNA_RNA_PACKAGE_ENERGY_PAR_H__
#define __VIENNA_RNA_PACKAGE_ENERGY_PAR_H__

#include "energy_const.h"

#define PUBLIC


extern double lxc37;   /* parameter for logarithmic loop
			  energy extrapolation            */

extern double stack37[NBPAIRS+1][NBPAIRS+1];
extern double stackdH[NBPAIRS+1][NBPAIRS+1]; /* stack enthalpies */
extern double entropies[NBPAIRS+1][NBPAIRS+1];  /* not used anymore */

extern double hairpin37[31];
extern double hairpindH[31];
extern double bulge37[31];
extern double bulgedH[31];
extern double internal_loop37[31];
extern double internal_loopdH[31];
extern double internal2_energy;
extern double old_mismatch_37[NBPAIRS+1][5][5];
extern double mismatchI37[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern double mismatchIdH[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern double mismatch1nI37[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern double mismatch23I37[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern double mismatch1nIdH[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern double mismatch23IdH[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern double mismatchH37[NBPAIRS+1][5][5];  /* same for hairpins */
extern double mismatchM37[NBPAIRS+1][5][5];  /* same for multiloops */
extern double mismatchHdH[NBPAIRS+1][5][5];  /* same for hairpins */
extern double mismatchMdH[NBPAIRS+1][5][5];  /* same for multiloops */
extern double mismatchExt37[NBPAIRS+1][5][5];
extern double mismatchExtdH[NBPAIRS+1][5][5];

extern double dangle5_37[NBPAIRS+1][5];      /* 5' dangle exterior of pair */
extern double dangle3_37[NBPAIRS+1][5];      /* 3' dangle */
extern double dangle3_dH[NBPAIRS+1][5];       /* corresponding enthalpies */
extern double dangle5_dH[NBPAIRS+1][5];

extern double int11_37[NBPAIRS+1][NBPAIRS+1][5][5]; /* 1x1 interior loops */
extern double int11_dH[NBPAIRS+1][NBPAIRS+1][5][5];

extern double int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5]; /* 2x1 interior loops */
extern double int21_dH[NBPAIRS+1][NBPAIRS+1][5][5][5];

extern double int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5]; /* 2x2 interior loops */
extern double int22_dH[NBPAIRS+1][NBPAIRS+1][5][5][5][5];

/* constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*(k-1) + ML_BASE*u  */
extern double ML_BASE37;
extern double ML_BASEdH;
extern double ML_closing37;
extern double ML_closingdH;
extern double ML_intern37;
extern double ML_interndH;

extern double TripleC37;
extern double TripleCdH;
extern double MultipleCA37;
extern double MultipleCAdH;
extern double MultipleCB37;
extern double MultipleCBdH;

/* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
/*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
extern double  MAX_NINIO;                   /* maximum correction */
extern double ninio37;
extern double niniodH;
/* penalty for helices terminated by AU (actually not GC) */
extern double TerminalAU37;
extern double TerminalAUdH;
/* penalty for forming bi-molecular duplex */
extern double DuplexInit37;
extern double DuplexInitdH;
/* stabilizing contribution due to special hairpins of size 4 (tetraloops) */
extern char Tetraloops[];  /* string containing the special tetraloops */
extern double  Tetraloop37[];  /* Bonus energy for special tetraloops */
extern double  TetraloopdH[];
extern char Triloops[];    /* string containing the special triloops */
extern double  Triloop37[]; /* Bonus energy for special Triloops */
extern double  TriloopdH[]; /* Bonus energy for special Triloops */
extern char Hexaloops[];    /* string containing the special triloops */
extern double  Hexaloop37[]; /* Bonus energy for special Triloops */
extern double  HexaloopdH[]; /* Bonus energy for special Triloops */

extern double GQuadAlpha37;
extern double GQuadAlphadH;
extern double GQuadBeta37;
extern double GQuadBetadH;

extern double Tmeasure;       /* temperature of param measurements */

#endif
