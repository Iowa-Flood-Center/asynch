#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>

int uncompress_gzfile(FILE *source, FILE *dest);
void zerr(int ret);

#endif

