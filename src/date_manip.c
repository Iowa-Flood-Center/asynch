#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>

#include "date_manip.h"


void copy_tm(struct tm *source,struct tm *dest)
{
	dest->tm_sec = source->tm_sec;
	dest->tm_min = source->tm_min;
	dest->tm_hour = source->tm_hour;
	dest->tm_mday = source->tm_mday;
	dest->tm_mon = source->tm_mon;
	dest->tm_year = source->tm_year;
	dest->tm_isdst = source->tm_isdst;
}

int days_in_month(int mon,int year)
{
	const int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	//int mon = current_time->tm_mon - 1;

	if(mon == 1)	//Damn you February...
	{
		//int year = current_time->tm_year;

		if(year % 4)	return days[mon];	//Not leap year
		if(year % 100)	return days[mon] + 1;	//Leap year
		if(year % 400)	return days[mon];	//Not leap year
		return days[mon] + 1;			//If this line is reached, then the current year is at least 2400.
							//I cannot explain to you in words how terrified I am to think that someone has
							//used code I wrote for 400 years. Seriously, something better and with
							//less bugs has surely come along. Also, anyone reading this in the year 2400 is
							//probably a robot, which is pretty cool.
							//--Scott,  April 2, 2014
	}
	else
		return days[mon];
}

