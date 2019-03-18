#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "LSRarefaction.h"

static int  verbose  = FALSE;

static int  nLines   = 10;

static char *usage[] = {"LSRarefactionS - \n",
                        "Required parameters:\n",
			"   -in   filename      parameter file\n",
			"   -c    float         desired coverage\n",
			"   -a    abundance file\n",
                        "Optional:\n",
			"   -b    integer       length of sample file to ignore\n",
			"   -s    integer       sampling frequency\n",
			"   -seed long          seed random number generator\n",
			"   -v                  verbose\n"};

double fMu(double x, void* pvParams)
{
  t_LSParams* ptLSParams = (t_LSParams *) pvParams;
  double dMDD            = ptLSParams->dMDash + x;
  double dLogP0          = logLikelihoodQuad(0, dMDD, ptLSParams->dV, ptLSParams->dNu);

  return (1.0 - exp(dLogP0)) - ptLSParams->dC;
}

double calcMu(double dC, t_LSParams *ptLSParams)
{
  double dLogMu = 0.0;

  ptLSParams->dC = dC;

  solveF(0, 1.0e7, fMu, ptLSParams, 1.0e-7, &dLogMu);

  return exp(dLogMu);
}

int compare_doubles(const void* a, const void* b) 
{
  double* arg1 = (double *) a;
  double* arg2 = (double *) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}   

int main(int argc, char* argv[]){
  t_Params    tParams;
  t_LSParams* atLSParams;
  int         i = 0, nSamples = 0;  
  double*     adMu = NULL;
  double dLower = 0.0, dMedian = 0.0, dUpper = 0.0;
  t_Data tData;

  gsl_set_error_handler_off();

  getCommandLineParams(&tParams, argc, argv);

  /*allocate memory for samples*/
  atLSParams = (t_LSParams *) malloc(MAX_SAMPLES*sizeof(t_LSParams));
  if(!atLSParams)
    goto memoryError;
  
  /*read in Monte-Carlo samples*/
  readSamples(&tParams, atLSParams, &nSamples);
  readAbundanceData(tParams.szAbundFile, &tData);
  
  adMu = (double *) malloc(sizeof(double)*nSamples);

  //printf("%d ",nSamples);
  for(i = 0; i < nSamples; i++){
    adMu[i] = ((double) tData.nJ)*calcMu(tParams.dCoverage, &atLSParams[i]);
    //printf("%f\n", adMu[i]);
    //fflush(stdout);
  }

  qsort(adMu, nSamples, sizeof(double), compare_doubles);


  dLower  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.025);
  dMedian = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.5);
  dUpper  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.975);

  printf("%.2e:%.2e:%.2e ", dLower, dMedian, dUpper);

  free(adMu);
  exit(EXIT_SUCCESS);
  
 memoryError:
  fprintf(stderr, "Failed to allocate memory in main aborting ...\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

int intCompare(const void *pvA, const void *pvB)
{
  int* pnA = (int *) pvA, *pnB = (int *) pvB;

  if(*pnA < *pnB)
    return +1;
  else if(*pnA == *pnB){
    return 0;
  }
  else{
    return -1;
  }
}

int doubleCompare(const void *pvA, const void *pvB)
{
  double* pnA = (double *) pvA, *pnB = (double *) pvB;

  if(*pnA < *pnB)
    return +1;
  else if(*pnA == *pnB){
    return 0;
  }
  else{
    return -1;
  }
}

void writeUsage(FILE* ofp)
{
  int i = 0;
  char *line;

  for(i = 0; i < nLines; i++){
    line = usage[i];
    fputs(line,ofp);
  }
}

char *extractParameter(int argc, char **argv, char *param,int when)
{
  int i = 0;

  while((i < argc) && (strcmp(param,argv[i]))){
    i++;
  }

  if(i < argc - 1){
    return(argv[i + 1]);
  }

  if((i == argc - 1) && (when == OPTION)){
    return "";
  }

  if(when == ALWAYS){
    fprintf(stdout,"Can't find asked option %s\n",param);
  }

  return (char *) NULL;
}

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[])
{
  char *szTemp = NULL;
  char *cError = NULL;

  /*get parameter file name*/
  ptParams->szInputFile  = extractParameter(argc,argv, INPUT_FILE,ALWAYS);  
  if(ptParams->szInputFile == NULL)
    goto error;

  ptParams->szAbundFile  = extractParameter(argc,argv, ABUND_FILE,ALWAYS);  
  if(ptParams->szAbundFile == NULL)
    goto error;

  /*get long seed*/
  szTemp  = extractParameter(argc,argv,SEED,OPTION);  
  if(szTemp != NULL){
    ptParams->lSeed = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->lSeed = 0;
  }

  /*get burn*/
  szTemp  = extractParameter(argc, argv, BURN, OPTION);  
  if(szTemp != NULL){
    ptParams->nBurn = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nBurn = DEF_BURN;
  }

  /*get coverage*/
  szTemp  = extractParameter(argc, argv, COVERAGE, ALWAYS);  
  if(szTemp != NULL){
    ptParams->dCoverage = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    goto error;
  }

 /*get long seed*/
  szTemp  = extractParameter(argc, argv, SAMPLE, OPTION);  
  if(szTemp != NULL){
    ptParams->nSample = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nSample = DEF_SAMPLE;
  }
 
  /*verbosity*/
  szTemp = extractParameter(argc, argv, VERBOSE, OPTION);
  if(szTemp != NULL){
    verbose = TRUE;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

void readSamples(t_Params *ptParams, t_LSParams *atLSParams, int *pnSamples)
{
  int  nSamples = 0;
  char *szInputFile = ptParams->szInputFile;
  char szLine[MAX_LINE_LENGTH];
  FILE *ifp = NULL;

  ifp = fopen(szInputFile, "r");

  if(ifp){
    while(fgets(szLine, MAX_LINE_LENGTH, ifp)){
      char *szTok = NULL, *szBrk = NULL, *pcError = NULL;
      int  nTime  = 0;

      /*remove trailing new line*/
      szBrk = strpbrk(szLine, "\n"); (*szBrk) = '\0';
      
      szTok = strtok(szLine, DELIM);

      nTime = strtol(szTok, &pcError, 10);
      if(*pcError != '\0') goto fileFormatError;

      if(nTime > ptParams->nBurn && nTime % ptParams->nSample == 0){
	
	szTok = strtok(NULL,DELIM);
	atLSParams[nSamples].dMDash = strtod(szTok, &pcError);
	if(*pcError != '\0') goto fileFormatError;

	szTok = strtok(NULL,DELIM);
	atLSParams[nSamples].dV = strtod(szTok, &pcError);
	if(*pcError != '\0') goto fileFormatError;

	szTok = strtok(NULL,DELIM);
	atLSParams[nSamples].dNu = strtod(szTok, &pcError);
	if(*pcError != '\0') goto fileFormatError;

	szTok = strtok(NULL,DELIM);
	atLSParams[nSamples].nS = strtol(szTok, &pcError, 10);
	if(*pcError != '\0') goto fileFormatError;

	nSamples++;
      }
    }

    fclose(ifp);
  }
  else{
    fprintf(stderr, "Failed to open file %s aborting\n", szInputFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  (*pnSamples) = nSamples;
  return;

 fileFormatError:
  fprintf(stderr, "Incorrectly formatted input file aborting\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

double f1(double x, void *pvParams)
{
  t_LSParams *ptLSParams = (t_LSParams *) pvParams;
  double dMDash = ptLSParams->dMDash, dV = ptLSParams->dV, dNu = ptLSParams->dNu;
  int n = ptLSParams->n;
  double t = ((x - dMDash)*(x - dMDash))/dV;
  double dExp  = x*((double) n) - exp(x);
  double dF  = pow(1.0 + t/dNu, -0.5*(dNu + 1.0));
  
  return exp(dExp)*dF;
}

double logStirlingsGamma(double dZ)
{
  return 0.5*log(2.0*M_PI) + (dZ - 0.5)*log(dZ) - dZ; 
}

double logLikelihoodQuad(int n, double dMDash, double dV, double dNu)
{
  gsl_integration_workspace *ptGSLWS = 
    gsl_integration_workspace_alloc(1000);
  double dLogFac1   = 0.0, dLogFacN  = 0.0;
  double dN = (double) n, dResult = 0.0, dError = 0.0, dPrecision = 0.0;
  gsl_function tGSLF;
  t_LSParams tLSParams;
  double dEst = dMDash + ((double)n)*dV, dA = 0.0, dB = 0.0;

  tLSParams.n = n; tLSParams.dMDash = dMDash; tLSParams.dV = dV; tLSParams.dNu = dNu;

  tGSLF.function = &f1;
  tGSLF.params   = (void *) &tLSParams;

  if(dNu < MAX_MU_GAMMA){
    dLogFac1 = gsl_sf_lngamma(0.5*(dNu + 1.0)) - gsl_sf_lngamma(0.5*dNu) - 0.5*log(M_PI*dNu);
  }
  else{
    dLogFac1 = 0.5*dNu*(log(0.5*(dNu + 1.0)) - log(0.5*dNu)) -0.5*log(2.0*M_PI) - 0.5;
  }
  
  if(n < 50){
    dLogFacN = gsl_sf_fact(n);
    dLogFacN = log(dLogFacN);
  }
  else if(n < 100){
    dLogFacN = gsl_sf_lngamma(dN + 1.0);
  }
  else{
    dLogFacN = logStirlingsGamma(dN + 1.0);
  }

  dA = -100.0; dB = 100.0;

  if(n < 10){
    dPrecision = HI_PRECISION;
  }
  else{
    dPrecision = LO_PRECISION;
  }

  gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
 
  //printf("%f %f\n", dResult, dError);

  gsl_integration_workspace_free(ptGSLWS);

  return log(dResult) - dLogFacN + dLogFac1 - 0.5*log(dV);
}

int solveF(double x_lo, double x_hi, double (*f)(double, void*), 
	   void* params, double tol, double *xsolve)
{
  int status, iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function F;
 
  F.function = f;
  F.params = params;
     
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
     
  do{
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
  
    status = gsl_root_test_interval (x_lo, x_hi, 0, tol);     
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  (*xsolve) = gsl_root_fsolver_root (s);
  gsl_root_fsolver_free (s);
     
  return status;
}


void readAbundanceData(const char *szFile, t_Data *ptData)
{
  int **aanAbund = NULL;
  int  i = 0, nNA = 0, nA = 0, nC = 0;
  int  nL = 0, nJ = 0;
  char szLine[MAX_LINE_LENGTH];
  FILE* ifp = NULL;

  ifp = fopen(szFile, "r");

  if(ifp){
    char* szTok   = NULL;
    char* pcError = NULL;

    fgets(szLine, MAX_LINE_LENGTH, ifp);

    szTok = strtok(szLine, DELIM2);
    
    nNA = strtol(szTok,&pcError,10);
    if(*pcError != '\0'){
      goto formatError;
    }
    
    aanAbund = (int **) malloc(nNA*sizeof(int*));

    for(i = 0; i < nNA; i++){
      aanAbund[i] = (int *) malloc(sizeof(int)*2);

      fgets(szLine, MAX_LINE_LENGTH, ifp);

      szTok = strtok(szLine, DELIM2);

      nA = strtol(szTok,&pcError,10);
      if(*pcError != '\0'){
	goto formatError;
      }

      szTok = strtok(NULL, DELIM2);

      nC = strtol(szTok,&pcError,10);
      if(*pcError != '\0'){
	goto formatError;
      }
      
      nL += nC;
      nJ += nC*nA;

      aanAbund[i][0]  = nA;
      aanAbund[i][1]  = nC;     
    }
  }
  else{
    fprintf(stderr, "Failed to open abundance data file %s aborting\n", szFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  ptData->nJ          = nJ;
  ptData->nL          = nL;
  ptData->aanAbund    = aanAbund;
  ptData->nNA         = nNA;
  return;

 formatError:
  fprintf(stderr, "Incorrectly formatted abundance data file\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}
