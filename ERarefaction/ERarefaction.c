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
#include "ERarefaction.h"

static char *usage[] = {"ERarefaction - Determines an emprical rarefraction curve by resampling\n",
                        "Required parameters:\n",
			"   -in   filename      parameter file  \n",
                        "Optional:\n",
			"   -sample int         how frequently to sample\n",
			"   -v                  verbose\n"};

static int  nLines   = 6;

static int  verbose  = FALSE;

double expectedDiversity(int n, t_Data *ptData)
{
  int i = 0;
  double dDenom = gsl_sf_lnchoose(ptData->nL, n);
  double dSum = 0.0;

  for(i = 0; i < ptData->nNA; i++){
    int nA = ptData->aanAbund[i][0];
    if(ptData->nL - nA >= n){

      double dNumer = gsl_sf_lnchoose(ptData->nL - nA, n);

      dSum += ((double) ptData->aanAbund[i][1])*exp(dNumer - dDenom);
    }

  }

  return ((double) ptData->nS) - dSum;
}

int main(int argc, char* argv[])
{
  int  i = 0, j = 0, r = 1;
  t_Params tParams;
  t_Data   tData;
  
  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);

  /*read in abundance distribution*/
  readAbundanceData(tParams.szInputFile, &tData);
     
  for(r = 1; r < tData.nL; r++){
    if(r % tParams.nSample == 0){
      printf("%d %f\n",r, expectedDiversity(r, &tData));
    }
  }
  
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
  szTemp  = extractParameter(argc,argv,SAMPLE,OPTION);  
  if(szTemp != NULL){
    ptParams->nSample = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nSample = 1;
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

void readAbundanceData(const char *szFile, t_Data *ptData)
{
  int **aanAbund = NULL;
  int  i = 0, nNA = 0, nA = 0, nC = 0;
  int  nL = 0, nS = 0;
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
      
      nS += nC;
      nL += nC*nA;

      aanAbund[i][0]  = nA;
      aanAbund[i][1]  = nC;     
    }
  }
  else{
    fprintf(stderr, "Failed to open abundance data file %s aborting\n", szFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  ptData->nS          = nS;
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

void freeAbundanceData(t_Data *ptData)
{
  int i = 0;

  for(i = 0; i < ptData->nNA; i++){
    free(ptData->aanAbund[i]);
  }
  free(ptData->aanAbund);
}
