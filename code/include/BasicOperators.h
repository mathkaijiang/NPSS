#ifndef __BasicOperators_h
#define __BasicOperators_h

#include "Head.h"


extern void printVecCplx(fftw_complex *src);
extern void printVecReal(double *src);
extern double normCplxInfty(fftw_complex *src, int n);
extern double normRealInfty(double *src, int n);
extern double minCplx(fftw_complex *src, int n);
extern void eliminate(fftw_complex *src, int n);
extern double mean_Cplx(fftw_complex *src, int n);

extern void printVect(double *src , int n);

extern double min_vector(double *src, int n);

extern double inpfp(fftw_complex *p, fftw_complex *q);
extern double nrmp(fftw_complex *p);

extern double inpfp(double *p, double *q);
extern double nrmp(double *p);

extern double rayleighq(fftw_complex *x, fftw_complex *v);

extern void FuncsLinear1Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1);

extern void FuncsLinear2Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2);

extern void FuncsLinear3Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3);

extern void FuncsLinear4Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4);

extern void FuncsLinear5Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4,
					  const double a5, const fftw_complex *F5);

extern void FuncAddToCplx(fftw_complex *rslt, int n,
			   const double a1, const fftw_complex *F1);

extern void FuncCplxAddAConst(fftw_complex *rslt, int n, const double a);
extern void FuncRealAddAConst(fftw_complex *rslt, int n, const double a);


extern void setCplxZero(fftw_complex rslt);
extern void setCplxZero_V(fftw_complex **W,int k,int cplxDofs);
extern void setCplxZero_v(fftw_complex *x,int cplxDofs);

extern void FuncsLinear1Real(double *rslt, int n,
					  const double a1, const double *F1);

extern void FuncsLinear2Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2);

extern void FuncsLinear3Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3);

extern void FuncsLinear4Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4);

extern void FuncAddToReal(double *rslt, int n,
				   const double a1, const double *F1);

extern void MatCopy(double **dst, double **src, int n, int m);
extern void MatsAdd(double **dst, double d, double **src, double s, int n, int m);
extern double MatsDot(double **lhs, double **rhs, int n, int m);
extern void MatPrint(double **matrix, int n, int m);
extern double innerProduct(double *lft, double *rht, int n);

extern void FFT_dot_constant(fftw_complex *rslt, int n,double a1);
extern void Dot_constant(double *rslt, int n,double a1);
extern void ngrad_cam(fftw_complex *phir, fftw_complex *gradr);
extern double ene_cam(fftw_complex *phi);
extern double ene_cam_F(fftw_complex *phi_F);
extern void hv_cam(fftw_complex *phir, fftw_complex *v, fftw_complex *hv);
extern void precond_camnew(fftw_complex *w);
extern void Orthoi_v2V(fftw_complex *v, fftw_complex **V,int k);
extern void Orthoi_Vi(fftw_complex **V,int k);


extern void get_Merge(double **UU,fftw_complex **V,fftw_complex **W,fftw_complex **Vp, int k); 
extern int get_nozero(int *emp,int k);
extern MatrixXd get_mat(double **UU,int *emp,int k,int m);
extern MatrixXd sort_eig(MatrixXd vector, MatrixXd value, int m);
extern void Printf_begin_end(fftw_complex *x);
extern void PV_Cplx(fftw_complex *v, fftw_complex **V,fftw_complex *tmpp,int k);
extern void hiosd_lobpcg_sieh_V0(fftw_complex *x, Option *option,fftw_complex **V0,int k0);

extern void copy_option_V(Option *option1,Option *option2,int k);
extern void copy_option_V(fftw_complex *v,Option *option2,int k);
extern void copy_option_V(fftw_complex *v,Option *option2,int k);
extern void copy_option_V(fftw_complex **V,Option *option2,int k);

#endif
