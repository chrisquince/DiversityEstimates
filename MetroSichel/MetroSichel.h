#ifndef LOG_NORMAL_H
#define LOG_NORMAL_H

#define DELIM " \n"
#define MAX_LINE_LENGTH   1024
#define MAX_WORD_LENGTH   128

#define SAMPLE_FILE_SUFFIX ".sample"

/*constants for simplex minimisation*/
#define PENALTY           1.0e20

#define INIT_A            0.1
#define INIT_B            1.0
#define INIT_G            -0.5
#define INIT_SIMPLEX_SIZE 0.1
#define MIN_SIMPLEX_SIZE  1.0e-5
#define MAX_SIMPLEX_ITER  100000

#define DEF_SIGMA     0.1
#define DEF_SIGMA_S   100.0

#define SLICE      10
#define PRECISION  1.0e-10

#define TRUE  1
#define FALSE 0


/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define INPUT_FILE       "-in"
#define OUT_FILE_STUB    "-out"
#define VERBOSE          "-v"
#define SAMPLE           "-s"
#define SIGMA_A          "-sigmaA"
#define SIGMA_B          "-sigmaB"
#define SIGMA_G          "-sigmaG"
#define SIGMA_S          "-sigmaS"
#define SEED             "-seed"

typedef struct s_Params
{
  long lSeed;

  char *szInputFile;

  char *szOutFileStub;

  double dSigmaA;
  
  double dSigmaB;

  double dSigmaG;

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

double logLikelihood(int n, double dAlpha, double dBeta, double dGamma);

int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params));

void freeAbundanceData(t_Data *ptData);

double negLogLikelihood(double dAlpha, double dBeta, double dGamma, int nS, void * params);

void outputResults(gsl_vector *ptX, t_Data *ptData);

void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams);

void* metropolis (void * pvInitMetro);

double Simpsons(double dA, double dB, int nN, double (*f1)(double x, void *pvParams), void* pvParams);

void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX);

#endif
