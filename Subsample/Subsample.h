#ifndef E_RAREFRACTION_H
#define E_RAREFRACTION_H

#define DELIM ",\n"
#define MAX_LINE_LENGTH   1024
#define MAX_WORD_LENGTH   128

#define TRUE  1
#define FALSE 0

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define INPUT_FILE       "-in"
#define OUT_FILE_STUB    "-out"
#define SAMPLE           "-i"
#define N_SAMPLES        "-r"
#define N_SIZE           "-n"
#define VERBOSE          "-v"
#define SEED             "-seed"

typedef struct s_Params
{
  long lSeed;

  int nI;

  int nR;

  int nN;	

  char *szInputFile;

  char *szOutFileStub;

} t_Params;




/*User defined functions*/

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readAbundanceData(int nI, const char *szFile, int **panN, int* pnS);

int compare_doubles(const void* a, const void* b);

#endif
