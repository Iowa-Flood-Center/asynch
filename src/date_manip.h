#ifndef DATE_MANIP_H
#define DATE_MANIP_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <time.h>

void copy_tm(struct tm *source,struct tm *dest);
int days_in_month(int mon,int year);

#endif

