#ifndef E_RAREFRACTION_H
#define E_RAREFRACTION_H

#define DELIM " \n"
#define MAX_LINE_LENGTH   1024
#define MAX_WORD_LENGTH   128

#define TRUE  1
#define FALSE 0

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define INPUT_FILE       "-in"
#define OUT_FILE_STUB    "-out"
#define VERBOSE          "-v"
#define SAMPLE           "-sample"

typedef struct s_Params
{
  long lSeed;

  int  nSample;

  char *szInputFile;

  char *szOutFileStub;

} t_Params;

typedef struct s_Data
{
  int nNA;
  
  int **aanAbund;

  int nS;

  int nL;
}t_Data;


/*User defined functions*/

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readAbundanceData(const char *szFile, t_Data *ptData);

int compare_doubles(const void* a, const void* b);

void freeAbundanceData(t_Data *ptData);

#endif
