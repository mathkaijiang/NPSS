#ifndef __Data_h
#define __Data_h

#include "Head.h"

typedef struct myVec
{
	double* Data;
	int* Index;
}mySortVec;

typedef struct myWave
{
	double* Wave;
	double Data;
}mySortWave;

extern char *model_initial_value;

extern int choice;
extern double angle;

extern double PARA_Q0, PARA_Q1, PARA_XI, PARA_lambda, PARA_tau, PARA_gamma;
extern int DIM, DimPhy, DimCpt, L;
extern int Nfold;
extern double **dirBox, **rcpBox;
extern int *NCpt;
extern int cplxDofs, realDofs;
extern double **ProjMatrix;
extern double *Gsquare;
extern double *epmckt2;
extern double *ikt2;
extern int **indKspace;
extern double **projPlane;
extern fftw_complex *rhoReal, *rhoCplx, *phi_F,*phi_R,*rhoReal1,*rhoCplx1;
extern fftw_complex *gradient, *gradientF;
extern fftw_complex *fftw_Ctmp, *fftw_Rtmp, *cplxTmp;
extern fftw_complex *quadTerm, *cubTerm, *quarTerm; 
extern fftw_complex *LaplaceR,*LaplaceF,*NonlinearR,*NonlinearF;
extern char *model_initial_value;
extern fftw_plan Planc2cFord, Planc2cBack,p,q;

struct Option
 {
     int k;
     int maxiter;
     int outputp;
     int outputd;     
     double dt;
     double epsf;
     double betat;
     double betau;

     fftw_complex **V;
 };
extern Option *option1,*option2;

#endif
