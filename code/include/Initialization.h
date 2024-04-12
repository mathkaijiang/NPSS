#ifndef __Initialization_h
#define __Initialization_h

#include "Head.h"

extern void initParameter(char system[20]);
extern void initialize();
extern void memAllocation();
extern void getProjPlane();
extern void getGsquare();
extern void initDenFourier();
extern void get_epmckt2();
extern void get_ikt2();
extern void get_initial_value(fftw_complex *r,int chice, double angle); 
#endif
