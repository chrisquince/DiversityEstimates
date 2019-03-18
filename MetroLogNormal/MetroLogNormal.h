#ifndef LOG_NORMAL_H
#define LOG_NORMAL_H

#define DELIM " \n"
#define MAX_LINE_LENGTH   1024
#define MAX_WORD_LENGTH   128

#define SAMPLE_FILE_SUFFIX ".sample"

/*constants for simplex minimisation*/
#define INIT_M_DASH       1.0
#define INIT_V            1.0
#define INIT_SIMPLEX_SIZE 1.0
#define INIT_S_SS         500.0
#define MIN_SIMPLEX_SIZE  1.0e-2

/*defaults for metropolis algorithm*/
#define DEF_SIGMA   0.1
#define DEF_SIGMA_S 100.0
#define SLICE       10

#define HI_PRECISION 1.0e-12
#define LO_PRECISION 1.0e-7

/*constants for calculated compound Poisson lognormal*/
#define V_MULT          25.0
#define MAX_QUAD        100
#define MAX_QUAD_DERIV  100

/*BOOLEANS*/
#define TRUE  1
#define FALSE 0

typedef struct s_LNParams
{
  double dMDash;

  double dV;

  int    n;
} t_LNParams;

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define INPUT_FILE       "-in"
#define OUT_FILE_STUB    "-out"
#define VERBOSE          "-v"
#define SAMPLE           "-s"
#define SIGMA_X          "-sigmaX"
#define SIGMA_Y          "-sigmaY"
#define SIGMA_S          "-sigmaS"
#define SEED             "-seed"

typedef struct s_Params
{
  long lSeed;

  char *szInputFile;

  char *szOutFileStub;

  double dSigmaX;
  
  double dSigmaY;

  double dSigmaS;
  
  int nIter;
} t_Params;

typedef struct s_Data
{
  int nNA;
  
  int **aanAbund;

  int nL;

  int nJ;
}t_Data;

typedef struct s_MetroInit
{
  t_Params *ptParams;

  t_Data   *ptData;

  gsl_vector* ptX;

  int nAccepted;

  long lSeed;

  int nThread;

} t_MetroInit;

/*User defined functions*/

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readAbundanceData(const char *szFile, t_Data *ptData);

double chao(t_Data *ptData);

int compare_doubles(const void* a, const void* b);

double f1(double x, void *pvParams);

double nLogLikelihood(const gsl_vector * x, void * params);

double logLikelihoodRampal(int n, double dMDash, double dV);

double logLikelihoodQuad(int n, double dMDash, double dV);

int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params));

void freeAbundanceData(t_Data *ptData);

double negLogLikelihood(double dMDash, double dV, int nS, void * params);

void outputResults(gsl_vector *ptX, t_Data *ptData);

void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, 
		 int* pnSDash, int nS, t_MetroInit *ptMetroInit);

void* metropolis (void * pvInitMetro);

double derivExponent(double x, void *pvParams);

int solveF(double ax, double bx, double (*f)(double, void*), void* params, double tol, double *xsolve);

void outputExpectations(FILE *ofp, double dMDash, double dV, int nS, t_Data *ptData);

void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX);

#endif
