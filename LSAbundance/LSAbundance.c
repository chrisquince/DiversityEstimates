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

#include "LSAbundance.h"

static int  verbose  = FALSE;

static int  nLines   = 9;

static char *usage[] = {"LSAbundance - \n",
                        "Required parameters:\n",
			"   -a    abundance data\n",
			"   -in   filename      parameter file  \n",
                        "Optional:\n",
			"   -b    integer       length of sample file to ignore\n",
			"   -s    integer       sampling frequency\n",
			"   -seed long          seed random number generator\n",
			"   -v                  verbose\n"};

void updateExpectations(double* adExpect, int nMax, double dMDash, double dV, double dNu, int nS, t_Data *ptData)
{
  int i = 0;

  for(i = 1; i <= nMax; i++){
    int nA = i;
    double dLog = 0.0, dP = 0.0;
    
    if(nA < MAX_QUAD){
      dLog = logLikelihoodQuad(nA, dMDash, dV, dNu);
    }
    else{
      dLog = logLikelihoodRampal(nA, dMDash, dV, dNu);
    }

    dP = exp(dLog);

    adExpect[i - 1]+= dP*nS;
  }
}


int main(int argc, char* argv[]){
  t_Params    tParams;
  t_LSParams* atLSParams;
  int         i = 0, nSamples = 0, nMax = 0;  
  t_Data      tData;
  double*     adExpect = NULL;

  gsl_set_error_handler_off();

  getCommandLineParams(&tParams, argc, argv);

  /*allocate memory for samples*/
  atLSParams = (t_LSParams *) malloc(MAX_SAMPLES*sizeof(t_LSParams));
  if(!atLSParams)
    goto memoryError;
  
  /*read in Monte-Carlo samples*/
  readSamples(&tParams, atLSParams, &nSamples);

  readAbundanceData(tParams.szAbundFile, &tData);

  nMax = tData.aanAbund[tData.nNA - 1][0];
  nMax = floor(pow(2.0,ceil(log((double) nMax)/log(2.0)) + 2.0) + 1.0e-7);

  adExpect = (double *) malloc(sizeof(double)*nMax);

  for(i = 0; i < nMax; i++){
    adExpect[i] = 0.0;
  }

  for(i = 0; i < nSamples; i++){
    updateExpectations(adExpect, nMax, 
		       atLSParams[i].dMDash, 
		       atLSParams[i].dV,
		       atLSParams[i].dNu,
		       atLSParams[i].nS, 
		       &tData);
  }

  for(i = 1; i <= nMax; i++){
    printf("%d %f\n",i,adExpect[i - 1]/((double) nSamples));
  }

  free(adExpect);
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

 /*get abundance file name*/
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

double f1Log(double x, void *pvParams)
{
  t_LNParams *ptLNParams = (t_LNParams *) pvParams;
  double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV;
  int n = ptLNParams->n;
  double dTemp = (x - dMDash);
  double dExp  = x*((double) n) - exp(x) - 0.5*((dTemp*dTemp)/dV);
  double dRet  = exp(dExp);

  return dRet;
}

double derivExponent(double x, void *pvParams)
{
  t_LNParams *ptLNParams = (t_LNParams *) pvParams;
  double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
  double dTemp = (x - dMDash)/dV, dRet = 0.0;

  dRet = ((double) n) - exp(x) - dTemp;

  return dRet;
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
  t_LNParams *ptLNParams = (t_LNParams *) params;
  F.function = f;
  F.params = params;
     
  //printf("%f %f %d %f %f\n",ptLNParams->dMDash, ptLNParams->dV, ptLNParams->n, x_lo, x_hi);
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

double logLikelihoodLNQuad(int n, double dMDash, double dV)
{
  gsl_integration_workspace *ptGSLWS = 
    gsl_integration_workspace_alloc(1000);
  double dLogFac1   = 0.0, dLogFacN  = 0.0;
  double dResult = 0.0, dError = 0.0, dPrecision = 0.0;
  gsl_function tGSLF;
  t_LNParams tLNParams;
  double dEst = dMDash + ((double)n)*dV, dA = 0.0, dB = 0.0;

  tLNParams.n = n; tLNParams.dMDash = dMDash; tLNParams.dV = dV;

  tGSLF.function = &f1Log;
  tGSLF.params   = (void *) &tLNParams;

  dLogFac1 = log(2.0*M_PI*dV);
  
  if(n < 50){
    dLogFacN = gsl_sf_fact(n);
    dLogFacN = log(dLogFacN);
  }
  else{
    dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
  }

  if(dEst > dV){
    double dMax = 0.0;
    double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
    double dVar   = 0.0;

    if(fabs(dUpper) > 1.0e-7){
      solveF(0.0, dUpper, derivExponent, (void *) &tLNParams, 1.0e-5, &dMax);
    }

    dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));

    dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
  }
  else{
    double dMax = 0.0;
    double dLower = dEst - dV;
    double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
    double dVar   = 0.0;

    if(fabs(dUpper - dLower) > 1.0e-7){
      solveF(dLower, dUpper, derivExponent, (void *) &tLNParams, 1.0e-5, &dMax);
    }
    else{
      dMax = 0.5*(dLower + dUpper);
    }
    dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));

    dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
  }
  
  if(n < 10){
    dPrecision = HI_PRECISION;
  }
  else{
    dPrecision = LO_PRECISION;
  }

  gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
 
  gsl_integration_workspace_free(ptGSLWS);

  return log(dResult) - dLogFacN -0.5*dLogFac1;
}

double logLikelihoodLNRampal(int n, double dMDash, double dV)
{
  double dN = (double) n;
  double dLogLik = 0.0, dTemp = gsl_pow_int(log(dN) - dMDash,2), dTemp3 = gsl_pow_int(log(dN) - dMDash,3);  

  dLogLik = -0.5*log(2.0*M_PI*dV) - log(dN) - (dTemp/(2.0*dV));

  dLogLik += log(1.0 + 1.0/(2.0*dN*dV)*(dTemp/dV + log(dN) - dMDash - 1.0) 
		 + 1.0/(6.0*dN*dN*dV*dV*dV)*(3.0*dV*dV - (3.0*dV - 2.0*dV*dV)*(dMDash - log(dN)) 
		 - 3.0*dV*dTemp + dTemp3));

  return dLogLik;
}

double logLikelihoodRampal(int n, double dMDash, double dV, double dNu)
{
  double dGamma = 0.5*(dNu + 1.0), dN = (double) n, dRN = 1.0/dN, dRSV = 1.0/(sqrt(dV)*sqrt(dNu));
  double dZ = (log(dN) - dMDash)*dRSV;
  double dDZDX = dRN*dRSV, dDZDX2 = -dRN*dRN*dRSV; 
  double dF = (1.0 + dZ*dZ);
  double dA = 0.0, dB = 0.0, dTemp = 0.0;
  double dLogFac1 = 0.0;

  if(dNu < MAX_MU_GAMMA){
    dLogFac1 = gsl_sf_lngamma(0.5*(dNu + 1.0)) - gsl_sf_lngamma(0.5*dNu) - 0.5*log(M_PI*dNu);
  }
  else{
    dLogFac1 = 0.5*dNu*(log(0.5*(dNu + 1.0)) - log(0.5*dNu)) -0.5*log(2.0*M_PI) - 0.5;
  } 
 
  dA = 4.0*dZ*dZ*dDZDX*dDZDX*dGamma*(dGamma + 1.0);
  dA /= dF*dF;

  dB = -2.0*dGamma*(dDZDX*dDZDX + dZ*dDZDX2);
  dB /= dF;

  dTemp = dRN + dA + dB;
 
  return -dGamma*log(dF) + log(dTemp) + dLogFac1 - 0.5*log(dV);
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
