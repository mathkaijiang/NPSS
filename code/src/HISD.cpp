#include "Data.h"
#include "Initialization.h"
#include "FftwToolkit.h"
#include "BasicOperators.h"
#include "DisplayResults.h"
#include "DualBox.h"
#include "LPsrc.h"
#include "HISD.h"
#include "Mytimer.h"
#include <Eigen/Dense>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace Eigen;
using namespace std;


void critical_point_HISD()
{
    int k = 4;
    hiosd_lobpcg_initialization(option2,k);
    option2->dt = 0.8;
    option2->betat = 8.0;
    option2->betau = 0.04;
    option2->epsf = 1e-7;
    option2->maxiter = 40000;
    option2->outputp = 100;
    option2->outputd = 1000;
    


    fftw_complex *v1,*v2,*v;
    v1   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	v2   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
    v   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
    get_initial_value(v1,2,89);
    get_initial_value(v2,2,-89);
    get_initial_value(v,2,90);
    for(int i=0;i<cplxDofs;i++) 
    {
        v[i][0] =  v[i][0] + v1[i][0] + v2[i][0];
        v[i][1] = 0.0;
    }


    printf("\n hhh1 \n");
    FuncsLinear2Cplx(rhoReal1, cplxDofs,1.0, rhoReal,0.15, v);
    //FuncsLinear2Cplx(rhoReal1, cplxDofs,1.0, rhoReal,0.05, option1->V[k-1]);    
    copy_option_V(option1,option2,k);
    printf("\n hhh2 \n");
   
    for(int i =0; i<cplxDofs;i++)
    {
        option2->V[k-1][i][0] = v[i][0];
    }
    
    hiosd_lobpcg_sieh(rhoReal1,option2);

    free(v);
    free(v1);
    free(v2);
}

void hiosd_lobpcg_sieh(fftw_complex *x, Option *option)
{
    int k,maxiter,outputd,outputp;
    k = option->k;
    maxiter = option->maxiter;
    outputp = option->outputp;
    outputd = option->outputd;

    double dt,epsf,betat,betau;
    dt = option->dt;
    epsf = option->epsf;
    betat = option->betat;
    betau = option->betau;

    double alpha[k];
    for(int i=0;i<k;i++) alpha[i] = 0.0;
    fftw_complex **U,**W,**Y,**Vp,**Up,**V;
    V = (fftw_complex **)malloc(sizeof(fftw_complex*) * k);
    U = (fftw_complex **)malloc(sizeof(fftw_complex*) * k);
    W = (fftw_complex **)malloc(sizeof(fftw_complex*) * k);
    Y = (fftw_complex **)malloc(sizeof(fftw_complex*) * k);
    Vp = (fftw_complex **)malloc(sizeof(fftw_complex*) * k);
    Up = (fftw_complex **)malloc(sizeof(fftw_complex*) * k);

    for(int i=0;i<k;i++)
    {
       V[i] = (fftw_complex *)malloc(sizeof(fftw_complex) *cplxDofs);
       U[i] = (fftw_complex *)malloc(sizeof(fftw_complex) *cplxDofs);
       W[i] = (fftw_complex *)malloc(sizeof(fftw_complex) *cplxDofs);
       Y[i] = (fftw_complex *)malloc(sizeof(fftw_complex) *cplxDofs);
       Vp[i]= (fftw_complex *)malloc(sizeof(fftw_complex) *cplxDofs);
       Up[i]= (fftw_complex *)malloc(sizeof(fftw_complex) *cplxDofs);
    }
    
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            V[i][j][0] = option->V[i][j][0];
            V[i][j][1] = 0.0;
        }
    }

    printf("\n rayleighq(x,v1) = %.3e ",rayleighq(rhoReal,V[0]) );
    printf("\n rayleighq(x,v2) = %.3e ",rayleighq(rhoReal,V[1]) );
    setCplxZero_V(W,k,cplxDofs);
    setCplxZero_V(Y,k,cplxDofs);
    setCplxZero_V(U,k,cplxDofs);
    setCplxZero_V(Vp,k,cplxDofs);
    setCplxZero_V(Up,k,cplxDofs);

    double **UU,**YY;
    UU = (double **)malloc(sizeof(double *) * 3*k);
    YY = (double **)malloc(sizeof(double *) * 3*k);
    for(int i=0;i<3*k;i++)
    {
        UU[i] = (double *)malloc(sizeof(double) *cplxDofs);
        YY[i] = (double *)malloc(sizeof(double) *cplxDofs);
    }


    int emp[3*k];
    int iter = 0,m;
    double nrmW ,temp1, nrmVpi;
    
    //正交化V
    for(int i=0;i<k;i++)
    {
        Orthoi_Vi(V,i);
        temp1 = nrmp(V[i]);
        FFT_dot_constant(V[i], cplxDofs,1.0/temp1);
    }
    
    fftw_complex *f,*g,*gp,*Dx,*tmpp,*tmppF,*nlnF,*nlnR,*xp,*Dg;
    f = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    gp = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    Dx = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    tmpp = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    tmppF = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    nlnR = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    nlnF = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    xp = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    g = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    Dg = (fftw_complex *)malloc(sizeof(fftw_complex)*cplxDofs);
    setCplxZero_v(tmpp,cplxDofs);
    setCplxZero_v(f,cplxDofs); 
    setCplxZero_v(Dx,cplxDofs);
    setCplxZero_v(tmppF,cplxDofs);
    setCplxZero_v(nlnR,cplxDofs);
    setCplxZero_v(nlnF,cplxDofs);
    setCplxZero_v(xp,cplxDofs);
    setCplxZero_v(g,cplxDofs);
    setCplxZero_v(Dg,cplxDofs);

    ngrad_cam(x, f);
    setCplxZero_v(gp,cplxDofs);
    PV_Cplx(f,V,tmpp,k);
    FuncsLinear2Cplx(Dx,cplxDofs,1.0,f,-2.0,tmpp);
    FuncRealAddAConst(Dx,cplxDofs,dt);

    double bta = 0.1,mean,res,maxf,minf;
    while(iter < maxiter)
    {
        for(int i=0;i<3*k;i++) emp[i] = 1;
        PV_Cplx(f,V,tmpp,k);
        FFT_dot_constant(tmpp,cplxDofs,2.0);//tmpp=2*PV(f)

        //printf("tmpp: %.3e ",nrmp(tmpp));
        //g = f-tmpp
        FuncsLinear2Cplx(g,cplxDofs,1.0,f,-1.0,tmpp);
        //Dg = g-gp
        FuncsLinear2Cplx(Dg,cplxDofs,1.0,g,-1.0,gp);
        //gp = g
        FuncsLinear1Cplx(gp,cplxDofs,1.0,g);
        bta = abs(inpfp(Dx,Dg)/inpfp(Dg,Dg));

        //printf("\n bta  %.3e  \n",bta );
        bta = min(bta,betat*dt);
        bta = max(bta,betau*dt);
        //printf("\n bta  %.3e  \n",bta );

        FuncsLinear1Cplx(xp,cplxDofs,1.0,x);//xp=x
        for(int i=0;i<cplxDofs;i++)
        {
            quadTerm[i][0] = x[i][0]*x[i][0];
            quadTerm[i][1] = 0.0;
            
            cubTerm[i][0] = quadTerm[i][0]*x[i][0];
            cubTerm[i][1] = 0.0;
        }

        p = fftw_plan_dft(DimCpt,NCpt,x,rhoCplx1,FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        FFT_dot_constant(rhoCplx1,cplxDofs,1.0/cplxDofs);
        setCplxZero(rhoCplx1[0]);

        //nln = PARA_tau +  PARA_gamma*quadTerm - cubTerm;
        for(int i=0;i<cplxDofs;i++)
        {
            nlnR[i][0] = PARA_tau*x[i][0] + PARA_gamma*quadTerm[i][0] - cubTerm[i][0];
            nlnR[i][1] = 0.0;
        }
        mean = mean_Cplx(nlnR, cplxDofs);
        //printf(" mean : %.3e",mean);
        FuncRealAddAConst(nlnR,cplxDofs,-mean);

        p = fftw_plan_dft(DimCpt,NCpt,nlnR,nlnF,FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        FFT_dot_constant(nlnF,cplxDofs,1.0/cplxDofs);

        p = fftw_plan_dft(DimCpt,NCpt,tmpp,tmppF,FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        FFT_dot_constant(tmppF,cplxDofs,1.0/cplxDofs);

        for(int i=0;i<cplxDofs;i++)
        {
            rhoCplx1[i][0]=( rhoCplx1[i][0]+bta*nlnF[i][0] - bta*tmppF[i][0])/(1+bta*epmckt2[i]);
            rhoCplx1[i][1]=0.0;
        }
        setCplxZero(rhoCplx1[0]);

        q = fftw_plan_dft(DimCpt,NCpt,rhoCplx1,x,FFTW_BACKWARD,FFTW_ESTIMATE);
        fftw_execute(q);
        
        //printf("\n nrmp(x) %.3e \n",nrmp(x));
        for(int i=0;i<cplxDofs;i++)
        {
            Dx[i][0] = x[i][0] - xp[i][0];
            Dx[i][1] = 0.0;
        }



        for(int i=0;i<k;i++)
        {           
            hv_cam(x, V[i], U[i]);
            alpha[i] = inpfp(V[i],U[i]);
            FuncsLinear2Cplx(W[i], cplxDofs,1, U[i],-alpha[i],V[i]);
            precond_camnew(W[i]);           
            Orthoi_v2V(W[i],V,k);
            Orthoi_Vi(W,i);
            nrmW = nrmp(W[i]);
            //printf("\n  nrmW[%d] = %.3e \n",i,nrmW);
            if(nrmW > 1e-8)
            {
                FFT_dot_constant(W[i], cplxDofs,1.0/nrmW);
                hv_cam(x, W[i], Y[i]);
            }
            else
            {
                setCplxZero_v(W[i],cplxDofs);
                setCplxZero_v(Y[i],cplxDofs);
                emp[k+i] = 0;
            }
        }

        for(int i=0;i<k;i++)
        {
            Orthoi_v2V(Vp[i],V,k);
            Orthoi_v2V(Vp[i],W,k);
            Orthoi_Vi(Vp,i);
            nrmVpi = nrmp(Vp[i]);
            //printf("\n  nrmVpi[%d] = %.3e \n",i,nrmVpi);
            
            if(nrmVpi > 1e-8)
            {
                FFT_dot_constant(Vp[i], cplxDofs,1.0/nrmVpi);
                hv_cam(x, Vp[i], Up[i]);
            }
            else
            {
                setCplxZero_v(Vp[i],cplxDofs);
                setCplxZero_v(Up[i],cplxDofs);
                emp[2*k+i] = 0;
            }
        }
        
        get_Merge(UU,V,W,Vp,k);
        get_Merge(YY,U,Y,Up,k);
        m = get_nozero(emp,k);
	    
        MatrixXd UUm = MatrixXd::Zero(m,cplxDofs);
        MatrixXd YYm = MatrixXd::Zero(m,cplxDofs);
        MatrixXd Pn = MatrixXd::Zero(m,m);
        
        UUm = get_mat(UU,emp,k,m);
        YYm = get_mat(YY,emp,k,m);
        Pn = UUm*YYm.transpose();
        Pn = (Pn.transpose() + Pn)*0.5;
        Pn = Pn/cplxDofs;        
        
        EigenSolver<MatrixXd> es(Pn);
        
        MatrixXd value= es.eigenvalues().real();
	    MatrixXd eta= es.eigenvectors().real();

        MatrixXd val_c = es.eigenvalues().imag();
        maxf = abs(  val_c.maxCoeff() );
        minf = abs(  val_c.minCoeff() );

        if(maxf > 1e-22 or minf > 1e-22 )
        {
            //printf("\n Complex characteristic value  \n" );
            break;
        }


        MatrixXd vector1 = eta;
        vector1 = sort_eig(eta,value,m);
        sort(value.data(),value.data()+value.size());
	    
        for(int i=0;i<k;i++)
        {
            for(int j=0;j<cplxDofs;j++)
            {
                Vp[i][j][0] = V[i][j][0];
                Vp[i][j][1] = 0.0;
                Up[i][j][0] = U[i][j][0];
                Up[i][j][1] = 0.0;
            }
        }
        
        MatrixXd eta1 = vector1.leftCols(k);
        MatrixXd Vm = eta1.transpose()*UUm;


        if(iter % 10 == 0)
        {
            cout << iter << " "<< value.transpose().leftCols(k)<<endl<<endl;

            //cout << "C  " << " "<< val_c.transpose()<<endl<<endl;
        }
        for(int i=0;i<k;i++)
        {
            for(int j=0;j<cplxDofs;j++)
            {
                V[i][j][0] = Vm(i,j);
                V[i][j][1] = 0.0;
            }
        }

	
        for(int i=0;i<k;i++)
        {
            temp1 = mean_Cplx(V[i], cplxDofs);
            for(int j=0;j<cplxDofs;j++)
            {
                V[i][j][0] = V[i][j][0] - temp1;
                V[i][j][1] = 0.0;
            }
        }

        for(int i=0;i<k;i++)
        {
            //Orthoi_Vi(V,i);
            Orthoi_v2V(V[i],V,i);
            temp1 = nrmp(V[i]);
            FFT_dot_constant(V[i], cplxDofs,1.0/temp1);
        }
 
        ngrad_cam(x,f);
        res = nrmp(f);
        
        if(iter % 10 == 0)
        {
            printf("\n %d \t bta = %.3e,  res = %.3e \n",iter,bta, res);
        }
        
        iter = iter +1;

        if(res < epsf)
        {
            break;
        }
    }

    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            option->V[i][j][0] = V[i][j][0];
            option->V[i][j][1] = 0.0;
        }
    }

    for(int i=0;i<k;i++)
    {
        fftw_free(V[i]);
        fftw_free(U[i]);
        fftw_free(W[i]);
        fftw_free(Y[i]);
        fftw_free(Up[i]);
        fftw_free(Vp[i]);
    }
    fftw_free(V);
    fftw_free(U);
    fftw_free(W);
    fftw_free(Y);
    fftw_free(Up);
    fftw_free(Vp);

    fftw_free(f);
    fftw_free(g);
    fftw_free(gp);
    fftw_free(Dx);
    fftw_free(tmpp);
    fftw_free(tmppF);
    fftw_free(xp);
    fftw_free(nlnF);
    fftw_free(nlnR);
    fftw_free(Dg);
}

