#include <cstdio> // NULL

#include "CStopWatch.h"

#ifdef WIN32

double CStopWatch::getElapsedTime() {	
    return (timer.stop - timer.start) / 1000. ;
}

#else

double CStopWatch::getElapsedTime() {	
	timeval res;
	timersub(&(timer.stop),&(timer.start),&res);
	return (res.tv_sec + res.tv_usec /1000000.0); // 10^6 uSec per second
}

#endif
