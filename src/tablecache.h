#ifndef TABLECACHE_H
#define TABLECACHE_H

#include <stddef.h>

int load_table(const char *fname, int version, void *table, size_t rank, size_t *dims, double *start, double *dx);
void write_table(const char *fname, int version, void *table, size_t rank, size_t *dims, double *start, double *dx);

#endif // TABLECACHE_H
