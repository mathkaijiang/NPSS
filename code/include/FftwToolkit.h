#ifndef __FftwToolkit_h
#define __FftwToolkit_h

#include "Head.h"

extern void getIndex(int **kspace, int dim);
extern int getIndex1D(int *k, int *Ndof, int dim);
extern void convolution(fftw_complex *src, fftw_complex *orig, int order);
extern void FftwC2C(fftw_complex *rslt, fftw_complex *src, char *direction);
extern void pad(fftw_complex *in, fftw_complex *out, int *Ndof, int *ExpandNdof, int dim);

#endif

