#define MODE_TURNER 0
#define MODE_UNCOUPLEDFTEMP 1
// #define MODE_UNIFORM 1
// #define MODE_NUSSINOV 2

#define DEFAULT_TEMP 37.0
#define DEBUG 0

typedef struct PartitionFunction {
   double **Z;
   double **ZB;
   double **ZM1;
   double **ZM;
} PartitionFunction;

PRIVATE short *encode_seq(const char *seq);
int structuralEntropyTurner(char* sequence, int verbose);
int structuralEntropyUncoupledFTemp(char* sequence, int verbose);
/* NUSSINOV AND UNIFORM NOT IMPLEMENTED 
int structuralEntropyUniform(char* sequence, int verbose);
int  structuralEntropyNussinov(char* sequence, int verbose);
int energyNJ(int i,int j, char* rna);
*/
