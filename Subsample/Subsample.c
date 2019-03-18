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

#include "Subsample.h"

static char *usage[] = {"Subsample - Subsamples an emprical rarefraction curve by resampling\n",
                        "Required parameters:\n",
			"   -in   filename      parameter file  \n",
			"   -i    integer       sample to subsample from\n",
			"   -r    integer       number of samples to generate\n",
			"   -n    integer       new sample size\n",
                        "Optional:\n",
			"   -seed long          seed random number generator\n",
			"   -v                  verbose\n"};

static int  nLines   = 9;

static int  verbose  = FALSE;

void setCP(double* adP, double *adC, int *anN, int nS, int nL)
{
  int i = 0;

  for(i = 0; i < nS; i++){
    adP[i] = ((double) anN[i])/((double) nL);
  }
  
  adC[i] = 0;
  for(i = 1; i <= nS; i++){
    adC[i] = adP[i - 1] + adC[i - 1];
  }
}

int sampleSpecies(gsl_rng *ptGSLRNG, double *adC, int nS)
{
  double dRand = gsl_rng_uniform(ptGSLRNG);
  int    nRet  = ((double) nS)*dRand;
  
  while(dRand <= adC[nRet]) nRet--;

  while(dRand > adC[nRet + 1]) nRet++;

  return nRet;
}

void writeSample(t_Params *ptParams, int nSample, int nR, int nS, int* anF)
{
  FILE *ofp = NULL;
  int i   = 0;
  char szFileName[MAX_LINE_LENGTH];
  
  sprintf(szFileName, "%s_%d_%d.dat",ptParams->szOutFileStub,nR, nSample);

  ofp = fopen(szFileName, "w");

  for(i = 0; i < nS; i++){
    fprintf(ofp,"%d %d\n",i, anF[i]);
  }
  
  fclose(ofp);
}

int main(int argc, char* argv[])
{
  int  i = 0, j = 0, r = 0;
  int  nCount = 0, nRS = 0, nD = 0;
  t_Params tParams;
  int      nS = 0, nL = 0;
  double   *adP = NULL;
  double   *adC = NULL;
  int      *anN = NULL;
  int      **aanF = NULL;
  gsl_rng            *ptGSLRNG = NULL;
  const gsl_rng_type *ptGSLRNGType = NULL;

  gsl_rng_env_setup();
     
  gsl_set_error_handler_off();
  
  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);

  ptGSLRNGType = gsl_rng_default;
  ptGSLRNG     = gsl_rng_alloc(ptGSLRNGType);

  gsl_rng_set(ptGSLRNG, tParams.lSeed);

  /*read in abundance distribution*/
  readAbundanceData(tParams.nI, tParams.szInputFile,&anN,&nS);

  for(i = 0; i < nS; i++){
    nL += anN[i];
  }

  adP = (double *) malloc(sizeof(double)*nS);
  adC = (double *) malloc(sizeof(double)*(nS + 1));
  aanF = (int    **) malloc(sizeof(int*)*nS);
  for(i = 0; i < nS; i++){
    aanF[i] = (int *) malloc(tParams.nR*sizeof(int));
  }

  setCP(adP, adC, anN, nS, nL);

  while(r < tParams.nR){
     
    for(i = 0; i < nS; i++){
      aanF[i][r] = 0;
    }

    for(i = 0; i < tParams.nN; i++){
      nD = sampleSpecies(ptGSLRNG, adC, nS);
      
      aanF[nD][r]++;
    }/*sample*/
	
    r++;
  }
  gsl_rng_free(ptGSLRNG);
	
  printf("OTU,");
  for(r = 0; r < tParams.nR - 1; r++){
    printf("r%d,",r);
  }
  printf("r%d\n",r);

  for(i = 0; i < nS; i++){
    printf("C%d,",i);
    for(r = 0; r < tParams.nR - 1; r++){
      printf("%d,",aanF[i][r]);
    }
    printf("%d\n",aanF[i][r]);
  } 

  for(i = 0; i < nS; i++){
    free(aanF[i]);
  }
  free(aanF);
  free(anN);
  free(adC);
  free(adP);

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
  //ptParams->szOutFileStub  = extractParameter(argc,argv,OUT_FILE_STUB,ALWAYS);  
  //if(ptParams->szOutFileStub == NULL)
    //goto error;

  szTemp  = extractParameter(argc,argv,SAMPLE,OPTION);  
  if(szTemp != NULL){
    ptParams->nI = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nI = 1;
  }

  szTemp  = extractParameter(argc,argv,N_SAMPLES,OPTION);  
  if(szTemp != NULL){
    ptParams->nR = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nR = 10;
  }

  szTemp  = extractParameter(argc,argv,N_SIZE,OPTION);  
  if(szTemp != NULL){
    ptParams->nN = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nN = 1000;
  }

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
  
  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

void readAbundanceData(int nI, const char *szFile, int **panN, int *pnS)
{
  int  *anN = NULL;
  int  i = 0, j = 0, nS = 0;
  char szLine[MAX_LINE_LENGTH];
  FILE* ifp = NULL;

  ifp = fopen(szFile, "r");

  if(ifp){
    char* szTok   = NULL;
    char* pcError = NULL;

    fgets(szLine, MAX_LINE_LENGTH, ifp);

    while(fgets(szLine, MAX_LINE_LENGTH, ifp) != NULL){
    	nS++;
    }
    fclose(ifp);

    ifp = fopen(szFile, "r");	
    fgets(szLine, MAX_LINE_LENGTH, ifp);	
    anN = (int *) malloc(nS*sizeof(int));

    for(i = 0; i < nS; i++){
    
      fgets(szLine, MAX_LINE_LENGTH, ifp);

      szTok = strtok(szLine, DELIM);

      for(j = 0; j < nI; j++){
 	szTok = strtok(NULL, DELIM);
      }
      anN[i] = strtol(szTok,&pcError,10);

      if(*pcError != '\0'){
	goto formatError;
      }
    }
  }
  else{
    fprintf(stderr, "Failed to open abundance data file %s aborting\n", szFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  (*pnS) = nS;
  (*panN) = anN;
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


