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
#include "MetroSichel.h"

static char *usage[] = {"MetroSichel - Fits the compound Poisson Sichel distn\n",
                        "Required parameters:\n",
			"   -out  filestub      output file stub\n",
			"   -in   filename      parameter file  \n",
                        "Optional:\n",
			"   -s    integer       generate integer MCMC samples\n",
			"   -seed   long        seed random number generator\n",
			"   -sigmaA float       std. dev. of alpha prop. distn\n",
			"   -sigmaB float       ...          beta             \n",
			"   -sigmaG float       ...          gamma            \n",
			"   -sigmaS float       ...          S                \n",
			"   -v                  verbose\n"};

static int  nLines   = 12;

static int  verbose  = FALSE;

int main(int argc, char* argv[])
{
  int  i = 0, nNA     = 0;
  t_Params tParams;
  t_Data   tData;
  gsl_vector* ptX = gsl_vector_alloc(4); /*parameter estimates*/

  
  gsl_rng_env_setup();
     
  gsl_set_error_handler_off();
  
  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);

  /*read in abundance distribution*/
  readAbundanceData(tParams.szInputFile, &tData);
     
  /*set initial estimates for parameters*/
  gsl_vector_set(ptX, 0, INIT_A);
  gsl_vector_set(ptX, 1, INIT_B);
  gsl_vector_set(ptX, 2, INIT_G);
  gsl_vector_set(ptX, 3, tData.nL*2);

  printf("D = %d L = %d Chao = %.2f\n",tData.nL, tData.nJ, chao(&tData));

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

  szTemp  = extractParameter(argc,argv,SIGMA_A,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaA = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaA = DEF_SIGMA;
  }
  
  szTemp  = extractParameter(argc,argv,SIGMA_B,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaB = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaB = DEF_SIGMA;
  }

  szTemp  = extractParameter(argc,argv,SIGMA_G,OPTION);  
  if(szTemp != NULL){
    ptParams->dSigmaG = strtod(szTemp,&cError);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dSigmaG = DEF_SIGMA;
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


double sd(int n, double dAlpha, double dBeta, double dGamma)
{
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

int bessel(double* pdResult, int n, double dAlpha, double dBeta, double dGamma)
{
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

double logLikelihood(int n, double dAlpha, double dBeta, double dGamma)
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

  status = bessel(&dRet,n, dAlpha,dBeta,dGamma);
  if(status == FALSE){
    dRet = sd(n, dAlpha,dBeta,dGamma);
  }
  
  return dRet - dLogFacN;
}

double nLogLikelihood(const gsl_vector * x, void * params)
{
  double dAlpha  = gsl_vector_get(x,0), dBeta = gsl_vector_get(x,1);
  double dGamma  = gsl_vector_get(x,2);
  int    nS = (int) floor(gsl_vector_get(x, 3));
  t_Data *ptData = (t_Data *) params;
  int    i       = 0;
  double dLogNot0 = 0.0, dLogL   = 0.0;
  double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
  
  if(dAlpha <= 0.0 || dBeta <= 0.0){
    return PENALTY;
  }

  for(i = 0; i < ptData->nNA; i++){
    double dLogP = 0.0;
    int    nA    = ptData->aanAbund[i][0];

    dLogP = logLikelihood(nA, dAlpha, dBeta, dGamma);
    
    dLogL += ((double) ptData->aanAbund[i][1])*dLogP;

    dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
    
  }

  dLog0 = logLikelihood(0, dAlpha, dBeta, dGamma);

  dLog1 = (nS - ptData->nL)*dLog0;

  dLog2 = - gsl_sf_lnfact(nS - ptData->nL);

  dLog3 = gsl_sf_lnfact(nS);
  
  dLogL += dLog1 + dLog2 + dLog3;

  /*return*/
  return -dLogL;
}

double negLogLikelihood(double dAlpha, double dBeta, double dGamma, int nS, void * params)
{
  t_Data *ptData = (t_Data *) params;
  int    i       = 0;
  double dLogNot0 = 0.0, dLogL   = 0.0;
  double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
  
  if(dAlpha <= 0.0 || dBeta <= 0.0){
    return PENALTY;
  }

  for(i = 0; i < ptData->nNA; i++){
    double dLogP = 0.0;
    int    nA    = ptData->aanAbund[i][0];

    dLogP = logLikelihood(nA, dAlpha, dBeta, dGamma);
    
    dLogL += ((double) ptData->aanAbund[i][1])*dLogP;

    dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
    
  }

  dLog0 = logLikelihood(0, dAlpha, dBeta, dGamma);

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
    gsl_vector_set(ss, i,INIT_SIMPLEX_SIZE*gsl_vector_get(ptX,i));
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
  double dDeltaA =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaA);
  double dDeltaB =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaB);
  double dDeltaG =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaG);
  int    nSDash = 0;

  gsl_vector_set(ptXDash, 0, gsl_vector_get(ptX,0) + dDeltaA);
  gsl_vector_set(ptXDash, 1, gsl_vector_get(ptX,1) + dDeltaB);
  gsl_vector_set(ptXDash, 2, gsl_vector_get(ptX,2) + dDeltaG);
  
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

  printf("\nML simplex: a = %.2f b = %.2f g = %.2f S = %.2f NLL = %.2f\n",dAlpha, dBeta, dGamma, dS, dL);
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

  printf("\nMCMC iter = %d sigmaA = %.2f sigmaB = %.2f sigmaG = %.2f sigmaS = %.2f\n",
	   ptParams->nIter, ptParams->dSigmaA, ptParams->dSigmaB, ptParams->dSigmaG, ptParams->dSigmaS);

  gsl_vector_memcpy(ptX1, ptX);

  gsl_vector_set(ptX2, 0, gsl_vector_get(ptX,0) + 2.0*ptParams->dSigmaA);
  gsl_vector_set(ptX2, 1, gsl_vector_get(ptX,1) + 2.0*ptParams->dSigmaB);
  gsl_vector_set(ptX2, 2, gsl_vector_get(ptX,2) + 2.0*ptParams->dSigmaG);   
  gsl_vector_set(ptX2, 3, gsl_vector_get(ptX,3) + 2.0*ptParams->dSigmaS); 

  gsl_vector_set(ptX3, 0, gsl_vector_get(ptX,0) - 2.0*ptParams->dSigmaA);
  gsl_vector_set(ptX3, 1, gsl_vector_get(ptX,1) - 2.0*ptParams->dSigmaB);
  gsl_vector_set(ptX3, 2, gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaG);
  if(gsl_vector_get(ptX,3) - 2.0*ptParams->dSigmaS > (double) ptData->nL){
        gsl_vector_set(ptX3, 3, gsl_vector_get(ptX,3) - 2.0*ptParams->dSigmaS);
  }
  else{
        gsl_vector_set(ptX3, 3, (double) ptData->nL);
  }
  
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
