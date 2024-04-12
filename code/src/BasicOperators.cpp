#include "Data.h"
#include "BasicOperators.h"
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
void printVecCplx(fftw_complex *src)
{
	for(int k = 0; k < cplxDofs; k++)
	{
		printf("src[%d][0] = %f \t src[%d][1] = %f\n", k, src[k][0], k, src[k][1]);
	}
}

void printVecReal(double *src)
{
	for(int k = 0; k < cplxDofs; k++)
	{
		printf("src[%d] = %f\n", k, src[k]);
	}
}

void printVect(double *src , int n)
{
	for(int k = 0; k < n; k++)
	{
		printf("src[%d] = %f\n", k, src[k]);
	}
}

double inpfp(fftw_complex *p, fftw_complex *q)
{
    double s = 0.0;
    for(int i=0;i<cplxDofs;i++)
    {
        s = s + p[i][0]*q[i][0];
    }
    s = s/cplxDofs;
    return s;
}

double inpfp(double *p, double *q)
{
    double s = 0.0;
    for(int i=0;i<cplxDofs;i++)
    {
        s = s + p[i]*q[i];
    }
    s = s/cplxDofs;
    return s;
}

double nrmp(fftw_complex *p)
{
    double s = 0.0;
    s = inpfp(p,p);
    s = sqrt(s);
    return s;
}

double nrmp(double *p)
{
    double s = 0.0;
    s = inpfp(p,p);
    s = sqrt(s);
    return s;
}

double rayleighq(fftw_complex *x, fftw_complex *v)
{
    double s1,s2,s;
    fftw_complex *Hv;
    Hv = (fftw_complex *)malloc( sizeof(fftw_complex) *cplxDofs );
    hv_cam(x,v,Hv);
    s1 = inpfp(v,Hv);
    s2 = inpfp(v,v);
    s = s1/s2;
    fftw_free(Hv);
    return s;
}

double normCplxInfty(fftw_complex *src, int n)
{
	double tmp;
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = src[i][0]*src[i][0]+src[i][1]*src[i][1];
		tmp = sqrt(tmp);
		rslt = (rslt > tmp ? rslt : tmp);
	}
	return rslt;
}

double normRealInfty(double *src, int n)
{
	double tmp;
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = src[i];
		rslt = (rslt > tmp ? rslt : tmp);
	}
	return rslt;
}

double mean_Cplx(fftw_complex *src, int n)
{
    double s =0;
    for(int i=0;i<n;i++)
    {
        s = s + src[i][0];
    }
    s = s/n;
    return s;
}

double minCplx(fftw_complex *src, int n)
{
	double tmp;
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = src[i][0]*src[i][0]+src[i][1]*src[i][1];
		tmp = sqrt(tmp);
		rslt = (rslt < tmp ? rslt : tmp);
	}
	return rslt;
}

double min_vector(double *src, int n)
{
	double tmp;
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp =  src[i];
		rslt = (rslt < tmp ? rslt : tmp);
	}
	return rslt;
}

void eliminate(fftw_complex *src, int n)
{
	for(int i = 0; i < n; i++)
	{
		if(fabs(src[i][0]) < 1.0e-15)
			src[i][0] = 0.0;
		if(fabs(src[i][1]) < 1.0e-15)
			src[i][1] = 0.0;
	}
}

void FuncsLinear1Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0];
		rslt[i][1] = a1*F1[i][1];
	}
}

void FuncsLinear2Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1];
	}
}

void FuncsLinear3Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0] + a3*F3[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1] + a3*F3[i][1];
	}
}

void FuncsLinear4Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0] + a3*F3[i][0] + a4*F4[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1] + a3*F3[i][1] + a4*F4[i][1];
	}
}

void FuncsLinear5Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4,
					  const double a5, const fftw_complex *F5)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0] + a3*F3[i][0] + a4*F4[i][0] + a5*F5[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1] + a3*F3[i][1] + a4*F4[i][1] + a5*F5[i][1];
	}
}

void FuncAddToCplx(fftw_complex *rslt, int n,
			   const double a1, const fftw_complex *F1)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] += a1*F1[i][0];
		rslt[i][1] += a1*F1[i][1];
	}
}

void FuncCplxAddAConst(fftw_complex *rslt, int n,
					   const double a)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] += a;
		rslt[i][1] += a;
	}
}

void FuncRealAddAConst(fftw_complex *rslt, int n,
					   const double a)
{
	for(int i = 0; i < n; i++)
		rslt[i][0] += a;
}


void setCplxZero(fftw_complex rslt)
{
	rslt[0] = 0.0; rslt[1] = 0.0;
}

void setCplxZero_V(fftw_complex **W,int k,int cplxDofs)
{
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            setCplxZero(W[i][j]);
        }
    }
}

void setCplxZero_v(fftw_complex *x,int cplxDofs)
{
    for(int i=0;i<cplxDofs;i++)
    {
        setCplxZero(x[i]);
    }
}


void FuncsLinear1Real(double *rslt, int n,
					  const double a1, const double *F1)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i];
}

void FuncsLinear2Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i];
}

void FuncsLinear3Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i] + a3*F3[i];
}

void FuncsLinear4Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i] + a3*F3[i] + a4*F4[i];
}

void FuncAddToReal(double *rslt, int n,
				   const double a1, const double *F1)
{
	for(int i = 0; i < n; i++)
		rslt[i] += a1*F1[i];
}

void MatCopy(double **dst, double **src, int n, int m)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			dst[i][j] = src[i][j];
}

void MatsAdd(double **dst, double d, double **src, double s, int n, int m)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			dst[i][j] = d*dst[i][j] + s*src[i][j];
}

double MatsDot(double **lhs, double **rhs, int n, int m)
{
	double rslt = 0;
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			rslt += lhs[i][j]*rhs[i][j];
	return rslt;
}

void MatPrint(double **matrix, int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			printf("matrix[%d][%d] = %.8f  ", i, j, matrix[i][j]);
		}
		printf("\n");
	}
}

double innerProduct(double *lft, double *rht, int n)
{
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
		rslt += lft[i]*rht[i];
	return rslt;
}

void FFT_dot_constant(fftw_complex *rslt, int n,
					   double a1)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*rslt[i][0];
		rslt[i][1] = a1*rslt[i][1];
	}
}

void Dot_constant(double *rslt, int n,double a1)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i] = a1*rslt[i];
	}
}



void ngrad_cam(fftw_complex *phir, fftw_complex *gradr) 
{
    fftw_complex *phif,*gradf;
    phif = (fftw_complex *) malloc(sizeof(fftw_complex) * cplxDofs);
    gradf= (fftw_complex *) malloc(sizeof(fftw_complex) * cplxDofs);
    
    p = fftw_plan_dft(DimCpt,NCpt,phir,phif,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(phif,cplxDofs,1.0/cplxDofs);
    setCplxZero(phif[0]);

    for(int i = 0; i < cplxDofs; i++)
	{   
        setCplxZero(gradf[i]);
		gradf[i][0] = epmckt2[i] * phif[i][0];
		gradf[i][1] = epmckt2[i] * phif[i][1];
	}
    
    q = fftw_plan_dft( DimCpt,NCpt,gradf,gradr,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(q);

    for(int i=0;i<cplxDofs;i++)
    {
        quadTerm[i][0] = phir[i][0]*phir[i][0];
        quadTerm[i][1] = 0.0;
        
        cubTerm[i][0] = quadTerm[i][0]*phir[i][0];
        cubTerm[i][1] = 0.0;
    }
    
    for(int i=0;i<cplxDofs;i++)
    {
        gradr[i][0] = gradr[i][0] - PARA_tau*phir[i][0] - \
                         PARA_gamma*quadTerm[i][0] + cubTerm[i][0];
        gradr[i][1] = 0.0;
    }

    //p = fftw_plan_dft(DimCpt,NCpt,gradr,gradf,FFTW_FORWARD,FFTW_ESTIMATE);
    //fftw_execute(p);
    //setCplxZero(gradf[0]);
    //FFT_dot_constant(gradf,cplxDofs,1.0/cplxDofs);

    //q = fftw_plan_dft(DimCpt,NCpt,gradf,gradr,FFTW_BACKWARD,FFTW_ESTIMATE);
    //fftw_execute(q);

    double mean = 0.0;
    mean = mean_Cplx(gradr, cplxDofs);
    FuncRealAddAConst(gradr,cplxDofs,-mean);

    FFT_dot_constant(gradr,cplxDofs,-1.0);

    fftw_free(phif);
    fftw_free(gradf);
}

double ene_cam(fftw_complex *phi)
{
    double E1,E2,E;
    p = fftw_plan_dft(DimCpt,NCpt,phi,phi_F,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(phi_F,cplxDofs,1.0/cplxDofs);
    setCplxZero(phi_F[0]);

    for(int i = 0; i < cplxDofs; i++)
	{   
        setCplxZero(gradientF[i]);
		fftw_Ctmp[i][0] = (pow(PARA_Q0,2)-Gsquare[i])*(pow(PARA_Q1, 2)-Gsquare[i])*phi_F[i][0];
		fftw_Ctmp[i][1] = (pow(PARA_Q0,2)-Gsquare[i])*(pow(PARA_Q1, 2)-Gsquare[i])*phi_F[i][1];
	}
    
    q = fftw_plan_dft( DimCpt,NCpt,fftw_Ctmp,fftw_Rtmp,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(q);

    for(int i=0;i<cplxDofs;i++)
    {
        LaplaceR[i][0] = 0.5*PARA_XI*fftw_Rtmp[i][0]*fftw_Rtmp[i][0];
        LaplaceR[i][1] = 0.0;
    }    
    p = fftw_plan_dft(DimCpt,NCpt,LaplaceR,LaplaceF,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(LaplaceF,cplxDofs,1.0/cplxDofs);
    E1 = LaplaceF[0][0];

    for(int i=0;i<cplxDofs;i++)
    {
        quadTerm[i][0] = phi[i][0]*phi[i][0];
        quadTerm[i][1] = 0.0;
        
        cubTerm[i][0] = quadTerm[i][0]*phi[i][0];
        cubTerm[i][1] = 0.0;

        quarTerm[i][0] = cubTerm[i][0]*phi[i][0];
        quarTerm[i][1] = 0.0;
    }
    
    for(int i=0;i<cplxDofs;i++)
    {
        NonlinearR[i][0] = -quadTerm[i][0]*PARA_tau/2-cubTerm[i][0]*PARA_gamma/3 + quarTerm[i][0]/4;
        NonlinearR[i][1] =  0.0;
    }

    p = fftw_plan_dft(DimCpt,NCpt,NonlinearR,NonlinearF,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(NonlinearF,cplxDofs,1.0/cplxDofs);

    E2 = NonlinearF[0][0];

    E = E1 + E2;
    return E;
}


double ene_cam_F(fftw_complex *phi_F)
{
    double E1,E2,E;
    setCplxZero(phi_F[0]);
    q = fftw_plan_dft(DimCpt,NCpt,phi_F,phi_R,FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(q);

    for(int i = 0; i < cplxDofs; i++)
	{   
        setCplxZero(fftw_Ctmp[i]);
		fftw_Ctmp[i][0] = (pow(PARA_Q0,2)-Gsquare[i])*(pow(PARA_Q1, 2)-Gsquare[i])*phi_F[i][0];
		fftw_Ctmp[i][1] = (pow(PARA_Q0,2)-Gsquare[i])*(pow(PARA_Q1, 2)-Gsquare[i])*phi_F[i][1];
	}
    
    q = fftw_plan_dft( DimCpt,NCpt,fftw_Ctmp,fftw_Rtmp,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(q);

    for(int i=0;i<cplxDofs;i++)
    {
        LaplaceR[i][0] = 0.5*PARA_XI*fftw_Rtmp[i][0]*fftw_Rtmp[i][0];
        LaplaceR[i][1] = 0.0;
    }    
    p = fftw_plan_dft(DimCpt,NCpt,LaplaceR,LaplaceF,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(LaplaceF,cplxDofs,1.0/cplxDofs);
    E1 = LaplaceF[0][0];
    for(int i=0;i<cplxDofs;i++)
    {
        quadTerm[i][0] = phi_R[i][0]*phi_R[i][0];
        quadTerm[i][1] = 0.0;
        
        cubTerm[i][0] = quadTerm[i][0]*phi_R[i][0];
        cubTerm[i][1] = 0.0;

        quarTerm[i][0] = cubTerm[i][0]*phi_R[i][0];
        quarTerm[i][1] = 0.0;
    }
    
    for(int i=0;i<cplxDofs;i++)
    {
        NonlinearR[i][0] = -quadTerm[i][0]*PARA_tau/2-cubTerm[i][0]*PARA_gamma/3 + quarTerm[i][0]/4;
        NonlinearR[i][1] =  0.0;
    }
    p = fftw_plan_dft(DimCpt,NCpt,NonlinearR,NonlinearF,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(NonlinearF,cplxDofs,1.0/cplxDofs);

    E2 = NonlinearF[0][0];
//    printf("\n E2 = %.5e",E2);
    E = E1 + E2;
    return E;
}

void hv_cam(fftw_complex *phir, fftw_complex *v, fftw_complex *hv)
{
    fftw_complex *phif,*phir2,*vf,*hvf;
    phif = (fftw_complex *) malloc( sizeof(fftw_complex) *cplxDofs );
    phir2= (fftw_complex *) malloc( sizeof(fftw_complex) *cplxDofs );
    hvf  = (fftw_complex *) malloc( sizeof(fftw_complex) *cplxDofs );
    vf   = (fftw_complex *) malloc( sizeof(fftw_complex) *cplxDofs );

    p = fftw_plan_dft(DimCpt,NCpt,v,vf,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(vf,cplxDofs,1.0/cplxDofs);

    for(int i=0;i<cplxDofs;i++)
    {
        hvf[i][0] = epmckt2[i]*vf[i][0];
        hvf[i][1] = epmckt2[i]*vf[i][1];
    }
    //setCplxZero(hvf[0]);
    q = fftw_plan_dft(DimCpt,NCpt,hvf,hv,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(q);
    
    for(int i=0;i<cplxDofs;i++)
    {
        hv[i][0] = hv[i][0] + (- PARA_tau - 2*PARA_gamma*phir[i][0] \
                   + 3*PARA_lambda*phir[i][0]*phir[i][0]) * v[i][0];
        hv[i][1] = 0.0;
    }
    
    double mean = 0.0;
    mean = mean_Cplx(hv, cplxDofs);
    //printf(" mean : %.3e",mean);
    FuncRealAddAConst(hv,cplxDofs,-mean);

    fftw_free(phif);
    fftw_free(phir2);
    fftw_free(vf);
    fftw_free(hvf);
}

void precond_camnew(fftw_complex *w)
{
    p = fftw_plan_dft(DimCpt,NCpt,w,fftw_Ctmp,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(fftw_Ctmp,cplxDofs,1.0/cplxDofs);
    for(int i=0;i<cplxDofs;i++)
    {
        fftw_Ctmp[i][0] = fftw_Ctmp[i][0] * ikt2[i];
        fftw_Ctmp[i][1] = fftw_Ctmp[i][1] * ikt2[i];
    }
    setCplxZero(fftw_Ctmp[0]);

    q = fftw_plan_dft(DimCpt,NCpt,fftw_Ctmp,w,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(q);
    double mean = 0.0;
    mean = mean_Cplx(w, cplxDofs);
    FuncRealAddAConst(w,cplxDofs,-mean);
}


void Orthoi_v2V(fftw_complex *v, fftw_complex **V,int k)
{
    double tmp = 0.0;
    for(int i=0;i<k;i++)
    {
        tmp = inpfp(v,V[i]);
        for(int j=0;j<cplxDofs;j++)
        {
            v[j][0] = v[j][0] -tmp*V[i][j][0];
            v[j][1] = 0.0;
        }
    }
}

void Orthoi_Vi(fftw_complex **V,int k)
{
    double tamp = 0.0;
    for(int i=0;i<k;i++)
    {
        tamp = inpfp(V[i],V[k]);
        for(int j=0;j<cplxDofs;j++)
        {
            V[k][j][0] = V[k][j][0]-tamp*V[i][j][0];
            V[k][j][1] = 0.0;
        }
    }
}

void get_Merge(double **UU,fftw_complex **V,fftw_complex **W,fftw_complex **Vp, int k)
{
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            UU[i][j] = V[i][j][0];
        }
    }

    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            UU[k+i][j] = W[i][j][0];
        }
    }

    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            UU[2*k+i][j] = Vp[i][j][0];
        }
    }
}

int get_nozero(int *emp,int k)
{
    int iter = 0;
    for(int i=0;i<3*k;i++)
    {
        if(emp[i]==1)
        {
            iter = iter + 1;
        }
    }
    return iter;
}

MatrixXd get_mat(double **UU,int *emp,int k,int m)
{
    MatrixXd UUm = MatrixXd::Zero(m,cplxDofs);
    int iter = 0;
    for(int i=0;i<3*k;i++)
    {
        if(emp[i] == 1)
        {
            for(int j=0;j<cplxDofs;j++)
                UUm(iter,j) = UU[i][j];
            iter = iter + 1;
        }
    }
    return UUm;
}

MatrixXd sort_eig(MatrixXd vector, MatrixXd value,int m)
{
    MatrixXd value1  = value;
    MatrixXd vector1 = vector;
    double maxi = value.maxCoeff();
    
    sort(value.data(),value.data()+value.size());
    for(int i=0;i<value.rows();i++)
    {
        for( int j=0;j<value.rows();j++)
        {
            if ( value(i) == value1(j) )
            {
                vector1.col(i) = vector.col(j);
                value1(j) = maxi+1;
                break;
            }           
        }
    }
    //cout << value << endl << endl;
    return vector1;
}


void Printf_begin_end(fftw_complex *x)
{
    printf("\n begin:  %.3e \t",x[0][0]);
    printf("end:  %.3e \n",x[cplxDofs-1][0]);
}

void PV_Cplx(fftw_complex *v, fftw_complex **V,fftw_complex *tmpp,int k)
{
    setCplxZero_v(tmpp,cplxDofs);
    double tmp;
    for(int i=0;i<k;i++)
    {
        tmp = inpfp(v,V[i]);
        for(int j=0;j<cplxDofs;j++)
        {
            tmpp[j][0] = tmpp[j][0] + tmp*V[i][j][0];
            tmpp[j][1] = 0.0;
        }
    }
}

void copy_option_V(Option *option1,Option *option2,int k)
{
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            option2->V[i][j][0] = option1->V[i][j][0];
            option2->V[i][j][1] = 0.0;
        }
    }
}

void copy_option_V(fftw_complex **V,Option *option2,int k)
{
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            option2->V[i][j][0] = V[i][j][0];
            option2->V[i][j][1] = 0.0;
        }
    }
}

void copy_option_V(fftw_complex *v,Option *option2,int k)
{
    for(int j=0;j<cplxDofs;j++)
    {
        option2->V[k][j][0] = v[j][0];
        option2->V[k][j][1] = 0.0;
    }
}




