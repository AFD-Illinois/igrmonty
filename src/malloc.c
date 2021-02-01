#include "decs.h"

void *malloc_rank1(int n1, int size)
{
  void *A;

  if((A = (void *)malloc(n1*size)) == NULL){
    fprintf(stderr,"malloc failure in malloc_rank1\n");
    exit(-1);
  }

  return A;
}

void **malloc_rank2(int n1, int n2, int size)
{
  void **A;

  if ((A = (void **)malloc(n1*sizeof(void *))) == NULL) {
    fprintf(stderr, "malloc_rank2 failure\n");
    exit(-1);
  }

  for (int i = 0; i < n1; i++) {
    A[i] = malloc_rank1(n2, size);
  }

  return A;
}

void ***malloc_rank3(int n1, int n2, int n3, int size)
{
  void ***A;

  if ((A = (void ***)malloc(n1*sizeof(void **))) == NULL) {
    fprintf(stderr, "malloc_rank3 failure\n");
    exit(-1);
  }

  for (int i = 0; i < n1; i++) {
    A[i] = malloc_rank2(n2, n3, size);
  }

  return A;
}

void ****malloc_rank4(int n1, int n2, int n3, int n4, int size)
{
  void ****A;

  if ((A = (void ****)malloc(n1*sizeof(void ***))) == NULL) {
    fprintf(stderr, "malloc_rank4 failure\n");
    exit(-1);
  }

  for (int i = 0; i < n1; i++) {
    A[i] = malloc_rank3(n2, n3, n4, size);
  }

  return A;
}

double **malloc_rank2_double(int n1, int n2)
{

  double **A;
  double *space;
  int i;

  space = malloc_rank1(n1*n2, sizeof(double));

  A = malloc_rank1(n1, sizeof(double *));

  if(!space || !A) {
    fprintf(stderr, "malloc_rank2_double failure\n");
    exit(-1);
  }

  for(i=0;i<n1;i++){
    A[i] = &(space[n2*i]);
    //A[i] = malloc_rank1(n2,sizeof(double *));
  }

  return A;
}

double ****malloc_rank4_double(int n1, int n2, int n3, int n4)
{

  double ****A;
  double *space;
  int i,j,k;

  space = malloc_rank1(n1*n2*n3*n4, sizeof(double));

  A = malloc_rank1(n1, sizeof(double *));

  if(!space || !A) {
    fprintf(stderr, "malloc_rank4_double failure\n");
    exit(-1);
  }

  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(double *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(double *));
      for(k=0;k<n3;k++){
        A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
      }
    }
  }

  return A;
}

void *****malloc_rank5(int n1, int n2, int n3, int n4, int n5, int size)
{
  void *****A;

  if ((A = (void *****)malloc(n1*sizeof(void ****))) == NULL) {
    fprintf(stderr, "malloc_rank5 failure\n");
    exit(-1);
  }

  for (int i = 0; i < N1; i++) {
    A[i] = malloc_rank4(n2, n3, n4, n5, size);
  }

  return A;
}
