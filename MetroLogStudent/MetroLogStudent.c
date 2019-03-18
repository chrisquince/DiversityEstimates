/*System includes*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/*GSL includes*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multimin.h>
#include <pthread.h>

/*User includes*/
#include "../Lib/FileUtils.h"
#include "../Lib/MatrixUtils.h"
#include "MetroLogStudent.h"

static char *usage[] = {"MetroLogStudent - Samples the Poisson Log-Student distn to SADs\n",
                        "Required parameters:\n",
			"   -out  filestub      output file stub\n",
			"   -in   filename      parameter file  \n",
                        "Optional:\n",
			"   -s    integer       generate integer mcmc samples\n",
			"   -seed   long        seed random number generator\n",
			"   -sigmaM float       std. dev. of mean prop. distn\n",
			"   -sigmaV float       \n",
			"   -sigmaN float       \n",
			"   -sigmaS float       \n",
			"   -v                  verbose\n"};

static int  nLines   = 11;

static int  verbose  = FALSE;

int main(int argc, char* argv[])
{
  int  i = 0, nNA     = 0;
  t_Params tParams;
  t_Data   tData;
  gsl_vector* ptX = gsl_vector_alloc(4); /*parameter estimates*/
  t_MetroInit atMetroInit[3];
  
  gsl_rng_env_setup();
     
  gsl_set_error_handler_off();
  
  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);

  /*read in abundance distribution*/
  readAbundanceData(tParams.szInputFile, &tData);
     
  /*set initial estimates for parameters*/
  gsl_vector_set(ptX, 0, INIT_M);
  gsl_vector_set(ptX, 1, INIT_V);
  gsl_vector_set(ptX, 2, INIT_N);
  gsl_vector_set(ptX, 3, tData.nL*2);

  printf("D = %d L = %d Chao = %f\n",tData.nL, tData.nJ, chao(&tData));

  minimiseSimplex(ptX, 4, (void*) &tData, &nLogLikelihood);

  outputResults(ptX, &tData);
   
  if(tParams.nIter > 0){
    mcmc(&tParams, &tData, ptX);
  }
  
  /*free up allocated memory*/
  gsl_vector_free(ptX);

  freeAbundanceData(&tData);

  exit(EXIT_SUCCESS);
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

  /*get out file stub*/
  ptParams->szOutFileStub  = extractParameter(argc,argv,OUT_FILE_STUB,ALWAYS);  
  if(ptParams->szOutFileStub == NULL)
    goto error;

  /*get out file stub*/
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
 
  /*verbosity*/
  szTemp = extractParameter(argc, argv, VERBOSE, OPTION);
  if(szTemp != NULL){
    verbose = TRUE;
  }

  szTemp  = extractParameter(argc,argv,SIGMA_M,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaM = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaM = DEF_SIGMA;
  }
  
  szTemp  = extractParameter(argc,argv,SIGMA_V,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaV = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaV = DEF_SIGMA;
  }

  szTemp  = extractParameter(argc,argv,SIGMA_N,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaN = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaN = DEF_SIGMA;
  }

  szTemp  = extractParameter(argc,argv,SIGMA_S,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaS = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaS = DEF_SIGMA_S;
  }

  szTemp  = extractParameter(argc,argv,SAMPLE,OPTION);  
  if(szTemp != NULL){
    ptParams->nIter = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nIter = 0;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
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

    szTok = strtok(szLine, DELIM);
    
    nNA = strtol(szTok,&pcError,10);
    if(*pcError != '\0'){
      goto formatError;
    }
    
    aanAbund = (int **) malloc(nNA*sizeof(int*));

    for(i = 0; i < nNA; i++){
      aanAbund[i] = (int *) malloc(sizeof(int)*2);

      fgets(szLine, MAX_LINE_LENGTH, ifp);

      szTok = strtok(szLine, DELIM);

      nA = strtol(szTok,&pcError,10);
      if(*pcError != '\0'){
	goto formatError;
      }

      szTok = strtok(NULL, DELIM);

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

int compare_doubles(const void* a, const void* b) 
{
  double* arg1 = (double *) a;
  double* arg2 = (double *) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}       

double chao(t_Data *ptData)
{
  double n1 = 0.0, n2 = 0.0;
  int **aanAbund = ptData->aanAbund;

  if(aanAbund[0][0] == 1 && aanAbund[1][0] == 2){
    n1 = (double) aanAbund[0][1]; n2 = (double) aanAbund[1][1];
  
    return ((double) ptData->nL) + 0.5*((n1*n1)/n2);
  }
  else{
    return -1.0;
  }
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


double nLogLikelihood(const gsl_vector * x, void * params)
{
  double dMDash  = gsl_vector_get(x,0), dV = gsl_vector_get(x,1);
  double dNu  = gsl_vector_get(x,2);
  int    nS = (int) floor(gsl_vector_get(x, 3));
  t_Data *ptData = (t_Data *) params;
  int    i       = 0;
  double dLogNot0 = 0.0, dLogL   = 0.0;
  double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
  
  if(dV <= 0.0 || dNu < 1.0){
    return PENALTY;
  }

  for(i = 0; i < ptData->nNA; i++){
    double dLogP = 0.0;
    int    nA    = ptData->aanAbund[i][0];

    if(nA < MAX_QUAD){
      dLogP = logLikelihoodQuad(nA, dMDash, dV, dNu);
    }
    else{
      dLogP = logLikelihoodRampal(nA, dMDash, dV, dNu);
    }
    
    dLogL += ((double) ptData->aanAbund[i][1])*dLogP;

    dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
    
  }

  dLog0 = logLikelihoodQuad(0, dMDash, dV, dNu);

  dLog1 = (nS - ptData->nL)*dLog0;

  dLog2 = - gsl_sf_lnfact(nS - ptData->nL);

  dLog3 = gsl_sf_lnfact(nS);
  
  dLogL += dLog1 + dLog2 + dLog3;

  /*return*/
  return -dLogL;
}

double negLogLikelihood(double dMDash, double dV, double dNu, int nS, void * params)
{
  t_Data *ptData = (t_Data *) params;
  int    i       = 0;
  double dLogNot0 = 0.0, dLogL   = 0.0;
  double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
  
  if(dV <= 0.0 || dNu < 1.0){
    return PENALTY;
  }

  for(i = 0; i < ptData->nNA; i++){
    double dLogP = 0.0;
    int    nA    = ptData->aanAbund[i][0];

    if(nA < MAX_QUAD){
      dLogP = logLikelihoodQuad(nA, dMDash, dV, dNu);
    }
    else{
      dLogP = logLikelihoodRampal(nA, dMDash, dV, dNu);
    }
    
    dLogL += ((double) ptData->aanAbund[i][1])*dLogP;

    dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
    
  }

  dLog0 = logLikelihoodQuad(0, dMDash, dV, dNu);

  dLog1 = (nS - ptData->nL)*dLog0;

  dLog2 = - gsl_sf_lnfact(nS - ptData->nL);

  dLog3 = gsl_sf_lnfact(nS);
  
  dLogL += dLog1 + dLog2 + dLog3;

  /*return*/
  return -dLogL;
}



int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params))
{
  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss;
  gsl_multimin_function minex_func;  
  size_t iter = 0;
  int i = 0, status;
  double size;

  /* Initial vertex size vector */
  ss = gsl_vector_alloc (nP);
     
  /* Set all step sizes to default constant */
  for(i = 0; i < nP; i++){
    gsl_vector_set(ss, i,INIT_SIMPLEX_SIZE*fabs(gsl_vector_get(ptX,i)));
  }

  /* Initialize method and iterate */
  minex_func.f = f;
  minex_func.n = nP;
  minex_func.params = pvData;
     
  s = gsl_multimin_fminimizer_alloc (T, nP);
  gsl_multimin_fminimizer_set(s, &minex_func, ptX, ss);
     
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
     
    if(status)
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, MIN_SIMPLEX_SIZE);
     
    if(status == GSL_SUCCESS){
      for(i = 0; i < nP; i++){
	gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
      }

      if(verbose) printf("converged to minimum at\n");
    }
    
    if(verbose){ 
      printf ("%5d ", iter);
    
      for (i = 0; i < nP; i++) printf("%10.3e ", gsl_vector_get(s->x, i));
    
      printf("f() = %7.3f size = %.3f\n", s->fval, size);
    }
  }
  while(status == GSL_CONTINUE && iter < MAX_SIMPLEX_ITER);
     
  for(i = 0; i < nP; i++){
    gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
  }

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

void freeAbundanceData(t_Data *ptData)
{
  int i = 0;

  for(i = 0; i < ptData->nNA; i++){
    free(ptData->aanAbund[i]);
  }
  free(ptData->aanAbund);
}

void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams)
{
  double dDeltaS =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaS);
  double dDeltaM =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaM);
  double dDeltaV =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaV);
  double dDeltaN =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaN);
  int    nSDash = 0;

  gsl_vector_set(ptXDash, 0, gsl_vector_get(ptX,0) + dDeltaM);
  gsl_vector_set(ptXDash, 1, gsl_vector_get(ptX,1) + dDeltaV);
  gsl_vector_set(ptXDash, 2, gsl_vector_get(ptX,2) + dDeltaN);
  
  //printf("%e %e %e\n",dDeltaA,dDeltaB,dDeltaG);

  nSDash = nS + (int) floor(dDeltaS);
  if(nSDash < 1){
    nSDash = 1;
  }
  (*pnSDash) = nSDash;
}


void outputResults(gsl_vector *ptX, t_Data *ptData)
{
  double dAlpha = 0.0, dBeta = 0.0, dGamma = 0.0, dS = 0.0, dL = 0.0;

  dAlpha = gsl_vector_get(ptX, 0);

  dBeta  = gsl_vector_get(ptX, 1);

  dGamma = gsl_vector_get(ptX, 2);
  
  dS = gsl_vector_get(ptX, 3);

  dL = nLogLikelihood(ptX, ptData);

  printf("\nML simplex: M = %.2f V = %.2f Nu = %.2f S = %.2f NLL = %.2f\n",dAlpha, dBeta, dGamma, dS, dL);
}

void* metropolis (void * pvInitMetro)
{
  t_MetroInit *ptMetroInit  = (t_MetroInit *) pvInitMetro;
  gsl_vector  *ptX          = ptMetroInit->ptX;
  t_Data      *ptData       = ptMetroInit->ptData;
  t_Params    *ptParams     = ptMetroInit->ptParams;
  gsl_vector  *ptXDash      = gsl_vector_alloc(4); /*proposal*/
  char *szSampleFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  const gsl_rng_type *T;
  gsl_rng            *ptGSLRNG;
  FILE    *sfp = NULL;
  int nS = 0, nSDash = 0, nIter = 0;
  double dRand = 0.0, dNLL = 0.0;
  void   *pvRet = NULL;

  /*set up random number generator*/
  T        = gsl_rng_default;
  ptGSLRNG = gsl_rng_alloc (T); 

  nS = (int) floor(gsl_vector_get(ptX,3));
  
  dNLL = negLogLikelihood(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), gsl_vector_get(ptX,2), nS,(void*) ptData);

  sprintf(szSampleFile,"%s_%d%s", ptParams->szOutFileStub, ptMetroInit->nThread, SAMPLE_FILE_SUFFIX);

  sfp = fopen(szSampleFile, "w");
  if(!sfp){
    exit(EXIT_FAILURE);
  }

  /*seed random number generator*/
  gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);

  /*now perform simple Metropolis algorithm*/
  while(nIter < ptParams->nIter){
    double dA = 0.0, dNLLDash = 0.0;

    getProposal(ptGSLRNG, ptXDash, ptX, &nSDash, nS, ptParams);
  
    dNLLDash = negLogLikelihood(gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), gsl_vector_get(ptXDash,2), nSDash, (void*) ptData);
    //printf("X' %e %e %e %d %f\n", gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), gsl_vector_get(ptXDash,2), nSDash, dNLLDash);
    //printf("X %e %e %e %d %f\n", gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), gsl_vector_get(ptX,2), nS, dNLL);
    dA = exp(dNLL - dNLLDash);
    if(dA > 1.0){
      dA = 1.0;
    }

    dRand = gsl_rng_uniform(ptGSLRNG);

    if(dRand < dA){
      gsl_vector_memcpy(ptX, ptXDash);
      nS = nSDash;
      dNLL = dNLLDash;
      ptMetroInit->nAccepted++;
    }
    
    if(nIter % SLICE == 0){
      fprintf(sfp, "%d,%.10e,%.10e,%.10e,%d,%f\n",nIter,gsl_vector_get(ptX, 0), gsl_vector_get(ptX, 1), gsl_vector_get(ptX, 2), nS, dNLL);    
      fflush(sfp);
    }

    nIter++;
  }

  fclose(sfp);

  /*free up allocated memory*/
  gsl_vector_free(ptXDash);
  free(szSampleFile);
  gsl_rng_free(ptGSLRNG);

  return pvRet;
}

void writeThread(t_MetroInit *ptMetroInit)
{
  gsl_vector *ptX = ptMetroInit->ptX;
    printf("%d: a = %.2f b = %.2f g = %.2f S = %.2f\n", ptMetroInit->nThread, 
	   gsl_vector_get(ptX, 0),
	   gsl_vector_get(ptX, 1),
	   gsl_vector_get(ptX, 2),
	   gsl_vector_get(ptX, 3));
}

void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX)
{
  pthread_t thread1, thread2, thread3;
  int       iret1  , iret2  , iret3;
  gsl_vector *ptX1 = gsl_vector_alloc(4), 
             *ptX2 = gsl_vector_alloc(4), 
             *ptX3 = gsl_vector_alloc(4);
  t_MetroInit atMetroInit[3];

  printf("\nMCMC iter = %d sigmaM = %.2f sigmaV = %.2f sigmaN = %.2f sigmaS = %.2f\n",
	   ptParams->nIter, ptParams->dSigmaM, ptParams->dSigmaV, ptParams->dSigmaN, ptParams->dSigmaS);

  gsl_vector_memcpy(ptX1, ptX);

  gsl_vector_set(ptX2, 0, gsl_vector_get(ptX,0) + 2.0*ptParams->dSigmaM);
  gsl_vector_set(ptX2, 1, gsl_vector_get(ptX,1) + 2.0*ptParams->dSigmaV);
  gsl_vector_set(ptX2, 2, gsl_vector_get(ptX,2) + 2.0*ptParams->dSigmaN);   
  gsl_vector_set(ptX2, 3, gsl_vector_get(ptX,3) + 2.0*ptParams->dSigmaS); 

  gsl_vector_set(ptX3, 0, gsl_vector_get(ptX,0) - 2.0*ptParams->dSigmaM);
  gsl_vector_set(ptX3, 1, gsl_vector_get(ptX,1) - 2.0*ptParams->dSigmaV);
  gsl_vector_set(ptX3, 2, gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaN);
  
  if(gsl_vector_get(ptX,3) - 2.0*ptParams->dSigmaS > (double) ptData->nL){
	gsl_vector_set(ptX3, 3, gsl_vector_get(ptX,3) - 2.0*ptParams->dSigmaS);   }
  else{
	gsl_vector_set(ptX3, 3, (double) ptData->nL);   }
  
  
  atMetroInit[0].ptParams = ptParams;
  atMetroInit[0].ptData   = ptData;
  atMetroInit[0].ptX      = ptX1;
  atMetroInit[0].nThread  = 0;
  atMetroInit[0].lSeed    = ptParams->lSeed;
  atMetroInit[0].nAccepted = 0;

  atMetroInit[1].ptParams = ptParams;
  atMetroInit[1].ptData   = ptData;
  atMetroInit[1].ptX      = ptX2;
  atMetroInit[1].nThread  = 1;
  atMetroInit[1].lSeed    = ptParams->lSeed + 1;
  atMetroInit[1].nAccepted = 0;

  atMetroInit[2].ptParams = ptParams;
  atMetroInit[2].ptData   = ptData;
  atMetroInit[2].ptX      = ptX3;
  atMetroInit[2].nThread  = 2;
  atMetroInit[2].lSeed    = ptParams->lSeed + 2;
  atMetroInit[2].nAccepted = 0;

  writeThread(&atMetroInit[0]);
  writeThread(&atMetroInit[1]);
  writeThread(&atMetroInit[2]);

  iret1 = pthread_create(&thread1, NULL, metropolis, (void*) &atMetroInit[0]);
  iret2 = pthread_create(&thread2, NULL, metropolis, (void*) &atMetroInit[1]);
  iret3 = pthread_create(&thread3, NULL, metropolis, (void*) &atMetroInit[2]);
  pthread_join(thread1, NULL);
  pthread_join(thread2, NULL);
  pthread_join(thread3, NULL);


  printf("%d: accept. ratio %d/%d = %f\n", atMetroInit[0].nThread, 
	 atMetroInit[0].nAccepted, ptParams->nIter,((double) atMetroInit[0].nAccepted)/((double) ptParams->nIter));

  printf("%d: accept. ratio %d/%d = %f\n", atMetroInit[1].nThread, 
	 atMetroInit[1].nAccepted, ptParams->nIter,((double) atMetroInit[1].nAccepted)/((double) ptParams->nIter));

  printf("%d: accept. ratio %d/%d = %f\n", atMetroInit[2].nThread,
	 atMetroInit[2].nAccepted, ptParams->nIter, ((double) atMetroInit[2].nAccepted)/((double) ptParams->nIter));
      
  gsl_vector_free(ptX1); gsl_vector_free(ptX2); gsl_vector_free(ptX3);
}
