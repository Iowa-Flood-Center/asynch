#ifndef COMPRESSION_H
#define COMPRESSION_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

int uncompress_gzfile(FILE *source, FILE *dest);
void zerr(int ret);

#endif

