#ifndef __DisplayResults_h
#define __DisplayResults_h

#include "Head.h"

extern void dispDensity(fftw_complex *rhoCplx, int step);
extern bool myComp(mySortVec a, mySortVec b);
extern bool scaleComp(mySortWave a, mySortWave b);
extern void dispCoeffFourier(fftw_complex *src, int step);
extern void dispPlaneWave(fftw_complex *src);
	
#endif
