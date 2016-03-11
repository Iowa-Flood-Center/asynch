#ifndef DATE_MANIP_H
#define DATE_MANIP_H

#include <stdio.h>
#include <time.h>

void copy_tm(struct tm *source,struct tm *dest);
int days_in_month(int mon,int year);

#endif

