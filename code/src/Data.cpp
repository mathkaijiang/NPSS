#include "Data.h"

int choice;
double angle;

double PARA_Q0, PARA_Q1, PARA_XI, PARA_lambda, PARA_tau, PARA_gamma;
int DIM, DimPhy, DimCpt;
int L;
int Nfold;
double **dirBox, **rcpBox;
int *NCpt;
int cplxDofs, realDofs;
double **ProjMatrix;
double *Gsquare;
double *epmckt2;
double *ikt2;
int **indKspace;
double **projPlane;

fftw_complex *rhoReal, *rhoCplx,*rhoReal1,*rhoCplx1;

fftw_complex *gradient;
fftw_complex *gradientF;
fftw_complex *fftw_Ctmp, *fftw_Rtmp, *cplxTmp;
fftw_complex *quadTerm, *cubTerm, *quarTerm;
fftw_complex *phi_F, *phi_R; 
fftw_complex *LaplaceR,*LaplaceF,*NonlinearR,*NonlinearF;
char *model_initial_value;
fftw_plan Planc2cFord, Planc2cBack,p,q;

Option *option1,*option2;
