#ifndef LOG_STUDENT_H
#define LOG_STUDENT_H

#define DELIM " \n"
#define MAX_LINE_LENGTH   1024
#define MAX_WORD_LENGTH   128

#define SAMPLE_FILE_SUFFIX ".sample"

/*constants for simplex minimisation*/
#define PENALTY           1.0e20

#define INIT_M            -10.0
#define INIT_V            20.0
#define INIT_N            20.0
#define INIT_SIMPLEX_SIZE 0.1
#define MIN_SIMPLEX_SIZE  1.0e-3
#define MAX_SIMPLEX_ITER  100000

#define DEF_SIGMA     0.1
#define DEF_SIGMA_S   100.0

#define SLICE      10
#define PRECISION  1.0e-10

#define TRUE  1
#define FALSE 0

typedef struct s_LSParams
{
  double dMDash;
  
  double dV;

  double dNu;

  int n;

} t_LSParams;


/*constants for calculated compound Poisson log-Student's*/
#define MAX_QUAD        100
#define HI_PRECISION    1.0e-12
#define LO_PRECISION    1.0e-7
#define MAX_MU_GAMMA    100

/*constants for calculated compound Poisson lognormal*/
#define V_MULT          25.0

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
#define SIGMA_M          "-sigmaM"
#define SIGMA_V          "-sigmaV"
#define SIGMA_N          "-sigmaN"
#define SIGMA_S          "-sigmaS"
#define SEED             "-seed"

typedef struct s_Params
{
  long lSeed;

  char *szInputFile;

  char *szOutFileStub;

  double dSigmaM;
  
  double dSigmaV;

  double dSigmaN;

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

double f1Log(double x, void *pvParams);

double nLogLikelihood(const gsl_vector * x, void * params);

double derivExponent(double x, void *pvParams);

double logLikelihoodQuad(int n, double dMDash, double dV, double dNu);

double logLikelihoodLNRampal(int n, double dMDash, double dV);

double logLikelihoodLNQuad(int n, double dMDash, double dV);

double logLikelihoodRampal(int n, double dMDash, double dV, double dNu);

int solveF(double x_lo, double x_hi, double (*f)(double, void*), 
	   void* params, double tol, double *xsolve);

int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params));

void freeAbundanceData(t_Data *ptData);

double negLogLikelihood(double dMDash, double dV, double dNu, int nS, void * params);

void outputResults(gsl_vector *ptX, t_Data *ptData);

void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams);

void* metropolis (void * pvInitMetro);

void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX);

#endif
