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

#include "IGRarefaction.h"

static int  verbose  = FALSE;

static int  nLines   = 10;

static char *usage[] = {"IGRarefaction - \n",
                        "Required parameters:\n",
			"   -in   filename      parameter file\n",
			"   -c    float         desired coverage\n",
			"   -a    abundance date file\n",
                        "Optional:\n",
			"   -b    integer       length of sample file to ignore\n",
			"   -s    integer       sampling frequency\n",
			"   -seed long          seed random number generator\n",
			"   -v                  verbose\n"};

double fMu(double x, void* pvParams)
{
  t_IGParams* ptIGParams = (t_IGParams *) pvParams;
  // double tx = x / 1667.0;
  
  double dAlphaDD = ptIGParams->dAlpha*sqrt(x);
  double dBetaDD  = ptIGParams->dBeta*x;
  double dLogP0   = logLikelihood(0, dAlphaDD, dBetaDD);

  return (1.0 - exp(dLogP0)) - ptIGParams->dC;
}

double calcMu(double dC, t_IGParams *ptIGParams)
{
  double dMu = 0.0;

  ptIGParams->dC = dC;

  solveF(1.0, 1.0e10, fMu, ptIGParams, 1.0e-7, &dMu);

  return dMu;
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
  t_IGParams* atIGParams;
  int         i = 0, nSamples = 0;  
  double*     adMu = NULL;
  double dLower = 0.0, dMedian = 0.0, dUpper = 0.0;
  t_Data tData;
  
  gsl_set_error_handler_off();

  getCommandLineParams(&tParams, argc, argv);

  /*allocate memory for samples*/
  atIGParams = (t_IGParams *) malloc(MAX_SAMPLES*sizeof(t_IGParams));
  if(!atIGParams)
    goto memoryError;
  
  /*read in Monte-Carlo samples*/
  readSamples(&tParams, atIGParams, &nSamples);

  readAbundanceData(tParams.szAbundFile, &tData);

  adMu = (double *) malloc(sizeof(double)*nSamples);

  //printf("%d ",nSamples);
  for(i = 0; i < nSamples; i++){
    adMu[i] = ((double) tData.nJ)*calcMu(tParams.dCoverage, &atIGParams[i]);
    //printf("%f\n", adMu[i]);
    //fflush(stdout);
  }

  qsort(adMu, nSamples, sizeof(double), compare_doubles);


  dLower  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.025);
  dMedian = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.5);
  dUpper  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.975);

  //printf("%.1f %.1f %.1f %.1f %.1f", dMedian, dUpper, dLower, dUpper, dLower);
  
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

  ptParams->szAbundFile  = extractParameter(argc,argv,ABUND_FILE,ALWAYS);  
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

void readSamples(t_Params *ptParams, t_IGParams *atIGParams, int *pnSamples)
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
	atIGParams[nSamples].dAlpha = strtod(szTok, &pcError);
	if(*pcError != '\0') goto fileFormatError;

	szTok = strtok(NULL,DELIM);
	atIGParams[nSamples].dBeta = strtod(szTok, &pcError);
	if(*pcError != '\0') goto fileFormatError;

	szTok = strtok(NULL,DELIM);
	atIGParams[nSamples].nS = strtol(szTok, &pcError, 10);
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

double fX(double x, double dA, double dB, double dNDash)
{
  double dTemp1 = (dA*(x - dB)*(x - dB))/x;

  return log(x) - (1.0/dNDash)*(x + dTemp1);
}

double f2X(double x, double dA, double dB, double dNDash)
{
  double dRet = 0.0, dTemp = 2.0*dA*dB*dB;

  dRet = (1.0/(x*x))*(1.0 + (1.0/dNDash)*(dTemp/x));

  return -dRet;
}


double sd(int n, double dAlpha, double dBeta)
{
  double dGamma = -0.5;
  double dA = 0.5*(-1.0 + sqrt(1.0 + (dAlpha*dAlpha)/(dBeta*dBeta)));
  double dN = (double) n, dNDash = dN + dGamma - 1.0, dRN = 1.0/dN;
  double dTemp1 = (0.5*dN)/(1.0 + dA), dTemp2 = 4.0*dRN*dRN*(1.0 + dA)*dA*dBeta*dBeta;
  double dXStar = dTemp1*(1.0 + sqrt(1.0 + dTemp2));
  double dFX = fX(dXStar, dA, dBeta, dNDash);
  double d2FX = -dNDash*f2X(dXStar, dA, dBeta, dNDash);
  double dLogK = 0.0, dGamma1 = dGamma;

  if(dGamma1 < 0.0){
    dGamma1 *= -1.0;
  }

  dLogK = gsl_sf_bessel_lnKnu(dGamma1,2.0*dA*dBeta);

  return -2.0*dA*dBeta -log(2.0) -dLogK -dGamma*log(dBeta) + dNDash*dFX + 0.5*log(2.0*M_PI) - 0.5*log(d2FX);
}

int bessel(double* pdResult, int n, double dAlpha, double dBeta)
{
  double dGamma  = -0.5;
  double dResult = 0.0;
  double dOmega = 0.0, dGamma2 = 0.0;
  double dLogK1 = 0.0, dLogK2 = 0.0;
  double dN = (double) n, dNu = dGamma + dN;
  double dTemp1 = 0.0;
 
  if(dNu < 0.0){
    dNu = -dNu;
  }

  if(dGamma < 0.0){
    dGamma2 = -dGamma;
  }
  else{
    dGamma2 = dGamma;
  }

  dOmega = sqrt(dBeta*dBeta + dAlpha*dAlpha) - dBeta;

  dLogK2 = gsl_sf_bessel_lnKnu(dNu, dAlpha);

  if(!gsl_finite(dLogK2)){
    if(dAlpha < 0.1*sqrt(dNu + 1.0)){
      //printf("l ");
      dLogK2 = gsl_sf_lngamma(dNu) + (dNu - 1.0)*log(2.0) - dNu*log(dAlpha);
    }
    else{
      //printf("sd ");
      (*pdResult) = dResult;
      return FALSE;
    }
  }
   
  dLogK1 = dGamma*log(dOmega/dAlpha) -gsl_sf_bessel_lnKnu(dGamma2,dOmega);
  
  dTemp1 = log((dBeta*dOmega)/dAlpha);

  dResult = dN*dTemp1 + dLogK2 + dLogK1;
  (*pdResult) = dResult;
  return TRUE;
}

double logLikelihood(int n, double dAlpha, double dBeta)
{
  double dLogFacN = 0.0;
  int status      = 0;
  double dRet     = 0.0;

  if(n < 50){
    dLogFacN = gsl_sf_fact(n);
    dLogFacN = log(dLogFacN);
  }
  else{
    dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
  }

  status = bessel(&dRet,n, dAlpha,dBeta);
  if(status == FALSE){
    dRet = sd(n, dAlpha,dBeta);
  }
  
  return dRet - dLogFacN;
}
