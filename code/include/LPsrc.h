#ifndef __LPsrc_h
#define __LPsrc_h

#include "Head.h"

extern void LPsystParameter();
extern double gradUpdate_LP(fftw_complex *phi,double tol);
extern double hamiltonUpdate_LP(double *timeSum);
extern void initialize_LP();
extern double gradientflow(fftw_complex *phir, double tol,int max_iter);
extern void hiosd_lobpcg_inip(fftw_complex *x ,Option *option1 );
extern void hiosd_lobpcg_initialization(Option *option1,int ki );
#endif
