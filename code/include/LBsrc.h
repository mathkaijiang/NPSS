#ifndef __LBsrc_h
#define __LBsrc_h

#include "Head.h"

extern void LBsystParameter();
extern double gradUpdate_LB(double tol);
extern double hamiltonUpdate_LB(double *timeSum);
extern void initialize_LB();

#endif
