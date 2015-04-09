/*---------------------------------------------
 * structuralEntropy.c
 *
 * J.A. Garcia-Martin - P.Clote
 * Program computes structural entropy for a given sequence S in 2 ways:
 * 1. Computing expected energy <E>(S) by dynamic programing:
 *       Q(S) = sum_over_all_structures(boltzman_factor(s) * E(s) 
 *       Z(S) = Partition function(S)
 *       <E>= Q(S)/Z(S)
 * 2. Computing expected energy <E>(S) estimating d/dT ln(Z(T)) for a given interval of T
 *       Q(S) = Partition function(S) uncoupling formal temperature and using T+delta_T as formal temperature.
 *       Z(S) = Partition function(S) 
 *       <E>= <E> = RT² * ((ln(Q(S))-ln(Z(S)))/delta_T)
 * Structural entropy H(s)=<E>/RT + ln(Z(S))
 * -----------------------------------------*/

// C headers
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h> //for INT_MAX
#include <stdlib.h>
#include <ctype.h>  //for toupper()
#include <string.h>


// Vienna headers
#include "fold.h"
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"


//Program headers
#include "convert_Vienna.h"
#include "RNAhairpin.h"
#include "myConst.h"
#include "misc.h"
#include "McCaskillSimple.h"
#include "structuralEntropy.h"


typedef long double DBL_TYPE;

void usage(char* name){
   printf ("\nUsage: %s \"sequence\" -s sequence -t temperature -d delta_Temp (compute structural entropy using <E> = RT² * d/dT ln(Z(T)) ) -z energy_is_zero [0|1] (default 0) -v (verbose, extended output)\n       %s -h (detailed help)  \n",name,name);
   printf ("Output:\n  - H     : Structural entropy\n  - <E>   : Expected energy\n  - Length: Sequence length\n  - G     : Ensemble free energy (-RT*ln(Z))\n");
   exit(1);
}


int main(int argc, char *argv[]){
	char sequence[MAXSIZE];
	sequence[0]='\0';
//	int mode = MODE_TURNER;
	int i;
	int verbose=0;
	double inputTemp = DEFAULT_TEMP;
	fTemp = DEFAULT_TEMP;
	double delta_T=0.0;

	/*---------------------------------------------------
	Check command line parameters
	---------------------------------------------------*/
	if (argc<2){ 
		usage(argv[0]);
		exit(1);
	}
	if (strcmp(argv[1],"-h")==0){ 
		printf("J.A. Garcia-Martin - P.Clote (2014)\n"
		       "Program computes structural entropy (H) for a given sequence S in 2 ways:\n"
		       "  1. Computing expected energy <E>(S) by dynamic programing:\n"
		       "      Q(S) = sum_over_all_structures(boltzman_factor(s) * E(s)\n"
		       "      Z(S) = Partition function(S)\n"
		       "      <E>  = Q(S)/Z(S)\n"
		       "  2. Computing expected energy <E>(S) estimating d/dT ln(Z(T)) for a given interval of T\n"
		       "      Q(S) = Partition function(S) uncoupling formal temperature and using T+delta_T as formal temperature.\n"
		       "      Z(S) = Partition function(S) \n"
		       "      <E>  =  RT² * ((ln(Q(S))-ln(Z(S)))/delta_T)\n"
		       "  Structural entropy H(s)=<E>/RT + ln(Z(S))\n\n"
		       "Input parameters are:\n"
		       "   <sequence>         : If FIRST argument is a valid nucleotide sequence it will be used input\n"
		       "   -s <sequence>      : Alternative flag for input sequence\n"
		       "   -t <temperature>   : Temperature in ºC\n"
		       "   -d <delta_T>       : Temperature variation used for estimating d/dT ln(Z(T))\n"
		       "                        NOTE: If this parameter is provided, H is computed estimating <E> = RT² * d/dT ln(Z(T))\n"
		       "   -v                 : Output includes method for computing H and the name of the ouput parameters\n"
		       "   -z <energy_is_zero>: [0|1] If value is 1, energies are set to 0, ouput is structural entropy for the uniform case \n"
		       "   -h                 : Print this message \n\n"
		       "Output format is:\n"
		       "   H: Structural entropy \t<E>: Expected energy\tLength: Sequence length\tG: Ensemble free energy (-RT*ln(Z))\n");
		usage(argv[0]);
	}
	int firstParameter=1;
	if (argv[1][0] != '-'){
		sequence[0]='@';
		strncpy(sequence+1,argv[1],strlen(argv[1])); 
		sequence[strlen(argv[1])+1]='\0'; //termination character for string
		firstParameter=2;
	}

	// Read optional parameters
	if (argc>2){
		for (i = firstParameter; i<argc; i++){
			if (argv[i][0] == '-'){
				switch (argv[i][1]){
					case 's': 
						if(++i<argc){
							sequence[0]='@';
							strncpy(sequence+1,argv[i],strlen(argv[i])); 
							sequence[strlen(argv[i])+1]='\0'; //termination character for string

							if ((sequence == NULL) || (sequence[0] == '-')){
								printf("\nError in sequence!\n");
								usage(argv[0]);
							}
						}
						else{
							printf("\nMissing input sequence after -s\n");
							usage(argv[0]);
						}
						break;
					case 'z': 
						if(++i<argc){
							if (sscanf(argv[i], "%d", &energy_is_zero)==0){
								printf("\nInvalid -z value !\n");
								usage(argv[0]);
							}
						}
						else{
							printf("\nMissing value after -z\n");
							usage(argv[0]);
						}

						break;
/* NUSSINOV AND UNIFORM NOT IMPLEMENTED 
					case 'm':
						if(++i<argc){
							if (sscanf(argv[i], "%d", &mode)==0)
								usage(argv[0]);
						}
						else{
							usage(argv[0]);
						}
						break;	
*/
					case 't':
						if(++i<argc){
							if (sscanf(argv[i], "%lf", &inputTemp)==0){
								printf("\nInvalid -t value !\n");
								usage(argv[0]);
							}
						}
						else{
							printf("\nMissing value after -t\n");
							usage(argv[0]);
						}
						break;
					case 'd':
						if(++i<argc){
							if (sscanf(argv[i], "%lf", &delta_T)==0){
								printf("\nInvalid -d value !\n");
								usage(argv[0]);
							}
						}
						else{
							printf("\nMissing value after -t\n");
							usage(argv[0]);
						}
						break;
					case 'v':
						verbose=1;
						break;
				}
			}
		}
	} 

	// Parameter validation
	if (strlen(sequence)==0){ 
		printf("No input sequence!\n");
		usage(argv[0]);
		exit(1);
	}
	else if(strlen(sequence)>MAXSIZE){
		printf("Maximum sequence length is %d.\n",MAXSIZE);
		exit(1);
	}
	if(inputTemp < -273){
		printf("Minimum temperature is -273.\n");
		exit(1);
	}
	if(energy_is_zero != 0 && energy_is_zero != 1){
		printf("Valid values for energy_is_zero are 0 and 1.\n");
		usage(argv[0]);
		exit(1);
	}
	if(delta_T < 0){
		printf("delta_T must be positive.\n");
		exit(1);
	}

	// Convert sequence to uppercase
	for(i=1;i<strlen(sequence);i++){
		sequence[i]=toupper(sequence[i]);
	}

	temperature=inputTemp;
	if(delta_T!=0){ // If delta T given compute structural entropy using <E> = RT² * d/dT ln(Z(T))
		fTemp = inputTemp+delta_T;
		return structuralEntropyUncoupledFTemp(sequence,verbose);
	}
	else{
		return structuralEntropyTurner(sequence,verbose);
	}


/* NUSSINOV AND UNIFORM NOT IMPLEMENTED 

	switch(mode){
		case MODE_UNIFORM:
			return structuralEntropyUniform(sequence,verbose);
			break;

		case MODE_NUSSINOV:
			return structuralEntropyNussinov(sequence,verbose);
			break;
		case MODE_UNCOUPLEDFTEMP:
			return structuralEntropyUncoupledFTemp(sequence,verbose);
			break;

		// Default mode
		case MODE_TURNER:
		default:
			return structuralEntropyTurner(sequence,verbose);
			break;
	}
*/
 
}

/*Below modified from fold.c*/
//Further modified by P.Clote on 12 May 2014
PRIVATE short *encode_seq(const char *seq) {
	unsigned int k,l;
	short *S0_out;
	l = strlen(seq);
	if ( (S0_out = (short *) calloc(1, sizeof(short)*(l+2) )) == NULL) {
		printf("Out of memory!\n");
		exit(1);
	}
	S0_out[0]= (short) l;
	for (k=1; k<=l; k++) { /* make numerical encoding of seq */
		S0_out[k]= (short) encode_char(toupper(seq[k-1]));
	}
	return S0_out;
}

int structuralEntropyTurner(char* sequence, int verbose){
	int i,j;
	// double totalpar; //total partition function
	// double totalNumStr; //total number of structures
	double expectedEnergy, expectedEnergyNorm,entropy;
 
	/*---------------------------------------------------
	Set up computation
	---------------------------------------------------*/

	CheckSequence(sequence);
	S0=encode_seq(sequence+1);
	seqlen=strlen(sequence)-1;
	Initialize_Params();
	make_pair_matrix();//needed for pair matching
	kT = (temperature+K0)*GASCONST/1000.0;
	fkT = (temperature+K0)*GASCONST/1000.0;

	//printf("%.15f\n",kT);
	if (energy_is_zero){
		ML_base=0;
		ML_close=0;
	}
	else{
		ML_base=(double)P->MLbase/100;
		ML_close=(double)P->MLclosing/100;
	}
  
	//allocate space for partition function 
	double **Z, **ZB, **ZM1, **ZM;
	Z   = Allocate2DMatrix( seqlen+1,seqlen+1);
	ZB  = Allocate2DMatrix( seqlen+1,seqlen+1);
	ZM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
	ZM  = Allocate2DMatrix( seqlen+1,seqlen+1);
	for(i=0;i<=seqlen;i++){ //necessary since Yang didn't use calloc
		for (j=0;j<=seqlen;j++) {
			if(j-i<4){
				ZB[i][j]=0; Z[i][j]=1; ZM1[i][j]=0; ZM[i][j]=0;
			}
		}
	}

	//allocate space for Q partition function 
	double **Q, **QB, **QM1, **QM;
	Q   = Allocate2DMatrix( seqlen+1,seqlen+1);
	QB  = Allocate2DMatrix( seqlen+1,seqlen+1);
	QM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
	QM  = Allocate2DMatrix( seqlen+1,seqlen+1);
	for(i=0;i<=seqlen;i++){ //necessary since Yang didn't use calloc
		for (j=0;j<=seqlen;j++) {
			QB[i][j]=0; Q[i][j]=0; QM1[i][j]=0; QM[i][j]=0;
		}
	}

	//compute partition function 
	mcCaskill(sequence, Z, ZB, ZM1, ZM);
	mcCaskillQ(sequence,Q,QB,QM1,QM, Z, ZB, ZM1, ZM);

	expectedEnergy = Q[1][seqlen]/Z[1][seqlen];
	expectedEnergyNorm = expectedEnergy/seqlen;
	entropy=(expectedEnergy/kT)+log(Z[1][seqlen]);
	double eEnsemble;
	eEnsemble=-kT*log(Z[1][seqlen]);

	if(verbose){
		printf("Computing expected energy <E>(S) by dynamic programing.\n");
		printf("H\t<E>\tLength\tG\n");
	}
	printf("%lf\t%lf\t%d\t%lf\n",entropy, expectedEnergy,seqlen,eEnsemble);
	return 0;
}

int structuralEntropyUncoupledFTemp(char* sequence, int verbose){
	int i,j;
	// double totalpar; //total partition function
	// double totalNumStr; //total number of structures
	double expectedEnergy, expectedEnergyNorm,entropy;
 
	/*---------------------------------------------------
	Set up computation
	---------------------------------------------------*/

	CheckSequence(sequence);
	S0=encode_seq(sequence+1);
	seqlen=strlen(sequence)-1;
	Initialize_Params();
	make_pair_matrix();//needed for pair matching
	kT = (temperature+K0)*GASCONST/1000.0;
	fkT = (fTemp+K0)*GASCONST/1000.0;

	//printf("%.15f\n",kT);
	if (energy_is_zero){
		ML_base=0;
		ML_close=0;
	}
	else{
		ML_base=(double)P->MLbase/100;
		ML_close=(double)P->MLclosing/100;
	}
  
	//allocate space for partition function 
	double **Z, **ZB, **ZM1, **ZM;
	Z   = Allocate2DMatrix( seqlen+1,seqlen+1);
	ZB  = Allocate2DMatrix( seqlen+1,seqlen+1);
	ZM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
	ZM  = Allocate2DMatrix( seqlen+1,seqlen+1);
	for(i=0;i<=seqlen;i++){ //necessary since Yang didn't use calloc
		for (j=0;j<=seqlen;j++) {
			if(j-i<4){
				ZB[i][j]=0; Z[i][j]=1; ZM1[i][j]=0; ZM[i][j]=0;
			}
		}
	}

	//allocate space for Q partition function using a different formal temperature
	double **Q, **QB, **QM1, **QM;
	Q   = Allocate2DMatrix( seqlen+1,seqlen+1);
	QB  = Allocate2DMatrix( seqlen+1,seqlen+1);
	QM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
	QM  = Allocate2DMatrix( seqlen+1,seqlen+1);
	for(i=0;i<=seqlen;i++){ //necessary since Yang didn't use calloc
		for (j=0;j<=seqlen;j++) {
			QB[i][j]=0; Q[i][j]=0; QM1[i][j]=0; QM[i][j]=0;
		}
	}

	fkT = (temperature+K0)*GASCONST/1000.0;
	//compute partition function 
	mcCaskill(sequence, Z, ZB, ZM1, ZM);

	fkT = (fTemp+K0)*GASCONST/1000.0;
	mcCaskill(sequence,Q,QB,QM1,QM);


	expectedEnergy = (GASCONST/1000.0)*(pow(temperature+K0,2))*((log(Q[1][seqlen])-log(Z[1][seqlen]))/(fTemp-temperature));

	expectedEnergyNorm = expectedEnergy/seqlen;
	entropy=(expectedEnergy/kT)+log(Z[1][seqlen]);
	double eEnsemble;
	eEnsemble=-kT*log(Z[1][seqlen]);

	if(verbose){
		printf("Computing expected energy by <E> = RT² * d/dT ln(Z(S,%fºC)) with delta_T=%e\n",temperature,(fTemp-temperature));		
		printf("H\t<E>\tLength\tG\n");
	}
	printf("%lf\t%lf\t%d\t%lf\n",entropy, expectedEnergy,seqlen,eEnsemble);
	return 0;
}

/* NUSSINOV AND UNIFORM NOT IMPLEMENTED 
int structuralEntropyUniform(char* sequence, int verbose){
}

int  structuralEntropyNussinov(char* sequence, int verbose){
}

int energyNJ(int i,int j, char* rna){
	//0<= i<j < n=len(rna)
	char ch1,ch2;

	if (DEBUG) return 0; //set energy to 0, so result same as uniform dist
		ch1 = toupper(rna[i]); ch2 = toupper(rna[j]);
	if ((ch1=='G' && ch2=='U') || (ch1=='U' && ch2=='G'))
		return(-1);
	//    return(0);
	if ((ch1=='A' && ch2=='U') || (ch1=='U' && ch2=='A'))
		return(-2);
	//    return(-1);
	//    return(0);
	if ((ch1=='G' && ch2=='C') || (ch1=='C' && ch2=='G'))
		return(-3);
	//    return(-1);
	//    return(0);
	return(INT_MAX);
}
*/
