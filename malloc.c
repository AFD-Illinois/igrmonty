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

  if ((A = (void **)malloc(n1*sizeof(void *))) == NULL) {
    fprintf(stderr, "malloc_rank2 failure\n");
    exit(-1);
  }

  for (int i = 0; i < N1; i++) {
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

  for (int i = 0; i < N1; i++) {
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

  for (int i = 0; i < N1; i++) {
    A[i] = malloc_rank3(n2, n3, n4, size);
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

