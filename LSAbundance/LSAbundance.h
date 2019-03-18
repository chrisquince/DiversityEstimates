#ifndef LOG_NORMAL_SAMPLE_H
#define LOG_NORMAL_SAMPLE_H


/*constants for calculated compound Poisson log-Student's*/
#define MAX_QUAD        100
#define HI_PRECISION    1.0e-12
#define LO_PRECISION    1.0e-7
#define MAX_MU_GAMMA    100

/*constants for calculated compound Poisson lognormal*/
#define V_MULT          25.0

typedef struct s_Params
{
  int nBurn;

  int nSample;

  long lSeed;

  char *szInputFile;

  char *szAbundFile;

  char *szOutFileStub;

} t_Params;

typedef struct s_LSParams
{
  int    nS;      /*number of species in community*/
  
  double dMDash;

  double dV;

  double dNu;

  int n;

} t_LSParams;

typedef struct s_Data
{
  int nNA;
  
  int **aanAbund;

  int nL;

  int nJ;
}t_Data;

typedef struct s_LNParams
{
  double dMDash;

  double dV;

  int    n;
} t_LNParams;


#define MAX_SAMPLES     1048576
#define MAX_LINE_LENGTH 1024
#define DELIM           ","
#define TRUE  1
#define FALSE 0
#define DELIM2 " \n"
/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define ABUND_FILE       "-a"
#define INPUT_FILE       "-in"
#define OUT_FILE_STUB    "-out"
#define BURN             "-b"
#define SAMPLE           "-s"
#define VERBOSE          "-v"
#define SEED             "-seed"

/*sampling parameters*/
#define DEF_BURN      2000000
#define DEF_SAMPLE    1000

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

int intCompare(const void *pvA, const void *pvB);

int doubleCompare(const void *pvA, const void *pvB);

void readSamples(t_Params *ptParams, t_LSParams *atLSParams, int *pnSamples);

double logLikelihoodRampal(int n, double dMDash, double dV, double dNu);

double logLikelihoodQuad(int n, double dMDash, double dV, double dNu);

void readAbundanceData(const char *szFile, t_Data *ptData);

void updateExpectations(double* adExpect, int nMax, double dMDash, double dV, double dNu, int nS, t_Data *ptData);

int solveF(double x_lo, double x_hi, double (*f)(double, void*), 
	   void* params, double tol, double *xsolve);

double derivExponent(double x, void *pvParams);

#endif
