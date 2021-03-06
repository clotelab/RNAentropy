#define MODE_TURNER 0
#define MODE_UNCOUPLEDFTEMP 1
// #define MODE_UNIFORM 1
// #define MODE_NUSSINOV 2

#define DEFAULT_TEMP 37.0

#define DEBUG 0

#define TURNER99 1999
#define TURNER04 2004
#define ANDRONESCU07 2007
#define DEFAULT_ENERGY_MODEL TURNER04

#define TURNER99_FILE "/param_files/rna_turner1999.par"
#define TURNER04_FILE "/param_files/rna_turner2004.par"
#define ANDRONESCU07_FILE "/param_files/rna_andronescu2007.par"

extern paramT *P;

typedef struct PartitionFunction {
   double **Z;
   double **ZB;
   double **ZM1;
   double **ZM;
} PartitionFunction;

PRIVATE short *encode_seq(const char *seq);
int structuralEntropyTurner(char* sequence, int verbose);
int structuralEntropyUncoupledFTemp(char* sequence, int verbose, int centered);
/* NUSSINOV AND UNIFORM NOT IMPLEMENTED 
int structuralEntropyUniform(char* sequence, int verbose);
int  structuralEntropyNussinov(char* sequence, int verbose);
int energyNJ(int i,int j, char* rna);
*/
