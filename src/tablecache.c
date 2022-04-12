#include "tablecache.h"
#include "hdf5_utils.h"

#include <stdlib.h>
#include <unistd.h>

// returns 0 on successful load
int load_table(const char *fname, int version, void *table, size_t rank, size_t *dims, double *start, double *dx)
{
  int rv = 0;
  if((rv = access(fname, F_OK)) != 0) {
    return rv;
  } 

  hdf5_open(fname);

  // check version is right
  int fversion = 0;
  hdf5_read_single_val(&fversion, "version", H5T_STD_I32LE);
  if (fversion != version) return -1;

  // check dimensions are right
  int frank = 0;
  hdf5_read_single_val(&frank, "rank", H5T_STD_I64LE);
  if (frank != rank) return -2;
  hsize_t rankcount[1] = { rank };
  hsize_t start1d[1] = { 0 };
  size_t *tempL = calloc(rank, sizeof(*tempL));
  hdf5_read_array(tempL, "dims", 1, rankcount, start1d, rankcount, rankcount, start1d, H5T_STD_I64LE);
  for (int i=0; i<rank; ++i) {
    if (dims[i] != tempL[i]) return -3;
  }
  free(tempL);

  // check start, dx correct
  double *tempD = calloc(rank, sizeof(*tempD));
  hdf5_read_array(tempD, "start", 1, rankcount, start1d, rankcount, rankcount, start1d, H5T_IEEE_F64LE);
  for (int i=0; i<rank; ++i) {
    if (start[i] != tempD[i]) return -4;
  }
  hdf5_read_array(tempD, "dx", 1, rankcount, start1d, rankcount, rankcount, start1d, H5T_IEEE_F64LE);
  for (int i=0; i<rank; ++i) {
    if (dx[i] != tempD[i]) return -5;
  }
  free(tempD);
  
  // load table
  hsize_t *start_zero = calloc(rank, sizeof(*start_zero));
  hsize_t *hdims = calloc(rank, sizeof(*hdims));
  for (int i=0; i<rank; ++i) hdims[i] = dims[i];
  hdf5_read_array(table, "table", rank, hdims, start_zero, hdims, hdims, start_zero, H5T_IEEE_F64LE);
  free(hdims);
  free(start_zero);

  hdf5_close();

  return 0;
}

void write_table(const char *fname, int version, void *table, size_t rank, size_t *dims, double *start, double *dx)
{
  hdf5_create(fname);

  // save version
  hdf5_write_single_val(&version, "version", H5T_STD_I32LE);

  // save dimensions
  hdf5_write_single_val(&rank, "rank", H5T_STD_I64LE);
  hsize_t rankcount[1] = { rank };
  hsize_t start1d[1] = { 0 };
  hdf5_write_array(dims, "dims", 1, rankcount, start1d, rankcount, rankcount, start1d, H5T_STD_I64LE);

  // save edges
  hdf5_write_array(start, "start", 1, rankcount, start1d, rankcount, rankcount, start1d, H5T_IEEE_F64LE);
  hdf5_write_array(dx, "dx", 1, rankcount, start1d, rankcount, rankcount, start1d, H5T_IEEE_F64LE);
 
  // save full table
  hsize_t *start_zero = calloc(rank, sizeof(*start_zero));
  hsize_t *hdims = calloc(rank, sizeof(*hdims));
  for (int i=0; i<rank; ++i) hdims[i] = dims[i];
  hdf5_write_array(table, "table", rank, hdims, start_zero, hdims, hdims, start_zero, H5T_IEEE_F64LE);
  free(hdims);
  free(start_zero);

  hdf5_close();
}

