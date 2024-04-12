#include "Data.h"
#include "Initialization.h"
#include "FftwToolkit.h"
#include "BasicOperators.h"
#include "DisplayResults.h"
#include "DualBox.h"
#include "LPsrc.h"
#include "Mytimer.h"
#include <Eigen/Dense>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace Eigen;
using namespace std;
void LPsystParameter()
{
// system parameters
    /*------------------ 1: qc, 2: C6, 3:CS6, 4: LQ , 5:LSQ , 6:T6 , 7: H , 8: LAm ----------*/
    choice = 1;  
    angle = 0.0;
    L = 112;
    DIM = 4;
	DimPhy = 2;
	DimCpt = 2;

	PARA_Q0 = 1.0;
	PARA_Q1 = 2.0*cos(PI/12.0);
	PARA_XI = 1;                 // penalty factor 
	PARA_tau    = 0.05;          // 二次项系数
	PARA_gamma  = 1.0;            // 三次项系数
	PARA_lambda = 1.0;            // 四次项系数

	Nfold = 12;  				  // 考虑结构的对称性
    

	printf("\t ***** System parameters: DIM = %d, Xi = %.3f, Tau = %.3f, Gamma = %.3f ***** \n", DimCpt, PARA_XI, PARA_tau, PARA_gamma);

	//model_initial_value="1fold";

	NCpt = (int *)malloc(sizeof(int)*DimCpt);
	for(int i = 0; i < DimCpt; i++) NCpt[i] = 1024;	

	//-------------------------------------------------------------------
	dirBox = (double **)malloc(sizeof(double*)*DimCpt);
	rcpBox = (double **)malloc(sizeof(double*)*DimCpt);
	for(int i = 0; i < DimCpt; i++)
	{
		dirBox[i] = (double *)malloc(sizeof(double)*DimCpt);
		rcpBox[i] = (double *)malloc(sizeof(double)*DimCpt);
	}
	for(int i = 0; i < DimCpt; i++)
	{
		for(int j = 0; j < DimCpt; j++)
		{
			dirBox[i][j] = 0.0;	rcpBox[i][j] = 0.0;
		}
	}
	double lenRecipBox = 1.0;
	for(int i = 0; i < DimCpt; i++)
    {
		rcpBox[i][i] = lenRecipBox;
        rcpBox[i][i] = rcpBox[i][i]/L;
    }
	getRecipLattice(rcpBox, dirBox, DimCpt);//???

	//-------------------------------------------------------------------
	ProjMatrix = (double **)malloc(sizeof(double*)*DimPhy);
	for(int i = 0; i < DimPhy; i++)
		ProjMatrix[i] = (double *)malloc(sizeof(double)*DIM);
	for(int i = 0; i < DimPhy; i ++)
		for(int j = 0; j < DIM; j ++)
			ProjMatrix[i][j] = 0.0;
	
//    for(int i = 0; i < DimPhy; i ++) ProjMatrix[i][i] = 1.0;
	for(int j = 0; j < DIM; j++)
	{
		ProjMatrix[0][j] = cos(2*PI*(j)/Nfold);
		ProjMatrix[1][j] = sin(2*PI*(j)/Nfold);
	}
	printf("\t\t\t\t\t === Projective Matrix === \n");
	MatPrint(ProjMatrix, DimPhy, DIM);
	printf("\n");
}

void initialize_LP()
{
	printf("\n\n \t\t\t\t ******* Initization ******* \n\n");
	mytimer_t timer;
// system parameters
	timer.reset();
	timer.start();
	LPsystParameter();
	timer.pause();
	printf("\t\t time cost of initialize parameters : %f seconds\n", timer.get_current_time());

// allocate memory 
	timer.reset();
	timer.start();
	memAllocation();
	timer.pause();
	printf("\t\t time cost of memory allocation : %f seconds\n", timer.get_current_time());

// 建立一维指标和物理空间指标的对应关系
	timer.reset();
	timer.start();
	getIndex(indKspace, DimCpt);
	timer.pause();
    printf(" \t\t getGsquare \n");
	printf("\t\t time cost of getIndex : %f seconds\n", timer.get_current_time());
    
// 平面波
	getProjPlane();
//  |G|^2
	getGsquare();

//  (q0-k^2)**2(q1-k^2)**2
    get_epmckt2();
    
//  预条件子
    get_ikt2();
//  密度rho的初值
	timer.reset();
	timer.start();
	//initDenFourier();
    get_initial_value(rhoReal,choice,angle);
	timer.pause();
	printf("\t\t time cost of initialize field : %f seconds\n", timer.get_current_time());
	
////  输出初始密度
//    dispCoeffFourier(rhoCplx, 0);
//    dispDensity(rhoCplx, 0);
}



double gradientflow(fftw_complex *phir, double tol,int max_iter)
{
    fftw_complex *phif,*phif1,*nlnF,*nlnR, *phir2,*phir3;
    phif = (fftw_complex *)malloc(sizeof(fftw_complex) * cplxDofs);
    nlnF = (fftw_complex *)malloc(sizeof(fftw_complex) * cplxDofs);
    nlnR = (fftw_complex *)malloc(sizeof(fftw_complex) * cplxDofs);
    phif1 = (fftw_complex *)malloc(sizeof(fftw_complex) * cplxDofs);
    phir2 = (fftw_complex *)malloc( sizeof(fftw_complex) * cplxDofs);
    phir3 = (fftw_complex *)malloc( sizeof(fftw_complex) * cplxDofs);

    p = fftw_plan_dft(DimCpt,NCpt,phir,phif,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(phif,cplxDofs,1.0/cplxDofs);
    setCplxZero(phif[0]);

    double res = 1,dt = 1,ene = 1,max_phi,max_nln;
    int iter = 0;

    while( res >  tol )
    {
        
         //计算非线性部分
        for(int i=0;i<cplxDofs;i++)
        {
            phir2[i][0] = phir[i][0]*phir[i][0];
            phir2[i][1] = 0.0;
    
            phir3[i][0]  = phir2[i][0]*phir[i][0];
            phir3[i][1]  = 0.0;
        }

        for(int i=0;i<cplxDofs;i++)
        {
            nlnR[i][0] = PARA_tau * phir[i][0] +  PARA_gamma * phir2[i][0] - PARA_lambda * phir3[i][0];
            nlnR[i][1] = 0.0;
        }
    
        p = fftw_plan_dft(DimCpt,NCpt,nlnR,nlnF,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_execute(p);
        FFT_dot_constant(nlnF,cplxDofs,1.0/cplxDofs);
    
        //迭代
        for(int i=0;i<cplxDofs;i++)
        {
            phif1[i][0] = (phif[i][0] + dt * nlnF[i][0])/(1+dt*epmckt2[i]);
            phif1[i][1] = (phif[i][1] + dt * nlnF[i][1])/(1+dt*epmckt2[i]);
        }        
        setCplxZero(phif1[0]);
        
        
        for(int i=0;i<cplxDofs;i++)
        {
            gradient[i][0] = (phif1[i][0] - phif[i][0])/dt;
            gradient[i][1] = (phif1[i][1] - phif[i][1])/dt;

            phif[i][0] = phif1[i][0];
            phif[i][1] = phif1[i][1];
        }

        q = fftw_plan_dft(DimCpt,NCpt,phif,phir,FFTW_BACKWARD,FFTW_ESTIMATE);
        fftw_execute(q);
        
        res = normCplxInfty(gradient, cplxDofs); 
        ene = ene_cam_F(phif);
        if(iter % 100 == 0)
            printf("\n iter = %d , res =  %.3e , ene = %.3e \n",iter,res,ene);
        //max_phi = normCplxInfty(phif,cplxDofs);
        //printf(" max_phi : %.3e",max_phi);

        iter = iter +1;

        if(iter > max_iter)
            break;
    }
    
    fftw_free(phif);
    fftw_free(phif1);
    fftw_free(nlnF);
    fftw_free(nlnR);
    fftw_free(phir2);
    fftw_free(phir3);
    return ene;
}


void hiosd_lobpcg_initialization(Option *option, int ki)
{
    option->k = ki;
    option->maxiter = 21;
    option1->dt = 0.1;
    option1->epsf = 1e-7; 
    option1->betat = 1.0;
    option1->betau = 1.0;
    option->outputp = 1;
    option1->outputd = 1;

    //printf("\n outputp:  %d",option1->outputp);
    option->V = (fftw_complex **) malloc( sizeof(fftw_complex*)*option->k );
    for(int i=0;i<option->k;i++)
    {
        option->V[i] = (fftw_complex*)malloc(sizeof(fftw_complex) *cplxDofs);
    }

    int k;
    k = option->k;

    MatrixXd Vm;
    Vm = MatrixXd::Random(k,cplxDofs);
    //Vm = MatrixXd::Identity(k,cplxDofs);

    double **V;
    V = (double **) malloc(sizeof(double *) *k);
    for(int i=0;i<k;i++)
        V[i] = (double *) malloc(sizeof(double ) *cplxDofs);
    

    for(int i=0;i<k;i++)
    {

        for(int j=0;j<cplxDofs;j++)
        {
            V[i][j] = Vm(i,j);
        }
    }
    
    double tmp,tmp1,s1;

    // 正交化 //
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<i;j++)
        {
            tmp = inpfp(V[i],V[j]);
            for(int l=0;l<cplxDofs;l++)
            {
                V[i][l] = V[i][l] - tmp*V[j][l];
            }
        }
        tmp1 = nrmp(V[i]);
        Dot_constant(V[i], cplxDofs,1.0/tmp1);
    }
    
    s1 = nrmp(V[k-1]);
    printf("\n  s1 :  %.3e \n",s1);    
    
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            option->V[i][j][0] = V[i][j];
            option->V[i][j][1] = 0.0;
        }
    }

    for(int i=0;i<k;i++)
    {
        free(V[i]);
    }
    free(V);
}

void hiosd_lobpcg_inip(fftw_complex *x ,Option *option )
{
    int k;
    k = option->k;
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
    int maxiter = option->maxiter;
    double nrmW ,temp1, nrmVpi;
    

    for(int i=0;i<k;i++)
    {
        Orthoi_Vi(V,i);
        temp1 = nrmp(V[i]);
        FFT_dot_constant(V[i], cplxDofs,1.0/temp1);
    }
 
   // printf("\n x: \n ");
   // Printf_begin_end(x);

   // printf("\n V[0]: \n");
   // Printf_begin_end(V[0]);

    double maxf=0,minf=0;
    while(iter < maxiter)
    {
        for(int i=0;i<3*k;i++) emp[i] = 1;
        for(int i=0;i<k;i++)
        {           
            hv_cam(x, V[i], U[i]);
            alpha[i] = inpfp(V[i],U[i]);
            //printf(" \n alpha[%d] = %.3e \n ",i,alpha[i]);
            FuncsLinear2Cplx(W[i], cplxDofs,1, U[i],-alpha[i],V[i]);
            //nrmW = nrmp(W[i]);
            //printf("\n  nrmW0[%d] = %.3e \n",i,nrmW);
 
            precond_camnew(W[i]);           

            //nrmW = nrmp(W[i]);
            //printf("\n  nrmW0[%d] = %.4e \n",i,nrmW);
 
            Orthoi_v2V(W[i],V,k);
            //nrmW = nrmp(W[i]);
            //printf("\n  nrmW0[%d] = %.3e \n",i,nrmW);
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
        //cout << "\n UU(0,0) : \n"<<UUm(0,0)<<endl<<endl;
        //cout<< "\n max UU  : \n"<<UUm.maxCoeff()<<endl<<endl;
        //cout<< "\n min UU  : \n"<<UUm.minCoeff()<<endl<<endl;
	    //cout<<"\n UU : \n" << UUm << endl<< endl;
        Pn = UUm*YYm.transpose();
        Pn = (Pn.transpose() + Pn)*0.5;
        Pn = Pn/cplxDofs;        
        //cout<<"\n Pn : \n" << Pn << endl<< endl;
        
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


        //cout << value<< endl << endl;

        //cout << " m = "<< m<< endl<< endl;
        MatrixXd vector1 = eta;
        vector1 = sort_eig(eta,value,m);
        sort(value.data(),value.data()+value.size());
	    //cout << value<< endl << endl;
	    
        //printf("\n eta : \n");
	    //cout << eta<< endl << endl;
	    
        //printf("\n vector1: \n ");
        //cout << vector1 << endl << endl;
        //eta = vector1;
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

        //printf("\n V[k-1]:  \n");
        //Printf_begin_end(V[k-1]);

//        printf("\n <V[0],V[k-1]> = %.3e",inpfp(V[0],V[k-1]) );
//        printf("\n <V[end],V[end]> = %.3e \n",inpfp(V[k-1],V[k-1]) );
	
	
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
        
       
        //printVect(alpha , k);
        //printf(" \n ");
        iter = iter +1;
    }
    for(int i=0;i<k;i++)
    {
        printf(" alpha[%d] = %.3e " ,i,rayleighq(x, V[i]) );
    }
    cout << endl;
 
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
}




