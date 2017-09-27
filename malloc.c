#include "decs.h"

void *malloc_rank1(int n1, int size)
{
  void *A;

  if((A = malloc(n1*size)) == NULL){
    fprintf(stderr,"malloc failure in malloc_rank1\n");
    exit(123);
  }

  return A;
}

void **malloc_rank2(int n1, int n2, int size)
{
  void **A;
  void *space;
  int i;

  space = malloc_rank1(n1*n2, size);

  A = malloc_rank1(n1, sizeof(void *));

  for(i = 0; i < n1; i++) A[i] = &(space[i*n2]);

  return A;
}

void ***malloc_rank3(int n1, int n2, int n3, int size)
{
  void ***A;
  void *space;
  int i,j;

  space = malloc_rank1(n1*n2*n3, size);

  A = malloc_rank1(n1, sizeof(void *));

  for(i = 0; i < n1; i++){
    A[i] = malloc_rank1(n2,sizeof(void *));
    for(j = 0; j < n2; j++){
      A[i][j] = &(space[n3*(j + n2*i)]);
    }
  }

  return A;
}

void ****malloc_rank4(int n1, int n2, int n3, int n4, int size)
{

  void ****A;
  void *space;
  int i,j,k;

  space = malloc_rank1(n1*n2*n3*n4, size);

  A = malloc_rank1(n1, sizeof(void *));

  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(void *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(void *));
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
  void *space;
  int i,j,k,l;

  space = malloc_rank1(n1*n2*n3*n4*n5, size);

  A = malloc_rank1(n1, sizeof(void *));

  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2, sizeof(void *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3, sizeof(void *));
      for(k=0;k<n3;k++){
        A[i][j][k] = malloc_rank1(n4, sizeof(void *));
        for(l=0;l<n4;l++){
          A[i][j][k][l] = &(space[n5*(l + n4*(k + n3*(j + n2*i)))]);
        }
      }
    }
  }

  return A;
}

