#include "Data.h"
#include "Initialization.h"
#include "FftwToolkit.h"
#include "BasicOperators.h"
#include "DisplayResults.h"
#include "Mytimer.h"
using namespace std;

void memAllocation()
{
	realDofs = cplxDofs = 1;
	for(int i = 0; i < DimCpt; i++)
	{
		realDofs *= NCpt[i];
		cplxDofs *= NCpt[i];
	}
	printf("\t\t\t\t\t === Discrete Modes === \n");
	for(int i = 0; i < DimCpt; i++)
	printf("\t NCpt[%d] = %d", i, NCpt[i]);
	printf("\n");
	printf("\t cplxDofs = %d,\t realDofs = %d\n\n", cplxDofs, realDofs);

	//-------------------------------------------------------------------
	printf("\t\t\t\t\t === Direction Box === \n");
	MatPrint(dirBox, DimCpt, DimCpt);
	printf("\n");
	printf("\t\t\t\t\t === reciprocal Box === \n");
	MatPrint(rcpBox, DimCpt, DimCpt);
	printf("\n");

	//-------------------------------------------------------------------
	indKspace = (int **)malloc(sizeof(int*)*cplxDofs);
	projPlane = (double **)malloc(sizeof(double*)*cplxDofs);
    	

	for(int i = 0; i < cplxDofs; i++)
	{
		indKspace[i] = (int *)malloc(sizeof(int)*DimCpt);
		projPlane[i] = (double *)malloc(sizeof(double)*DimPhy);
	}

    
	rhoCplx   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	rhoReal   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
    rhoCplx1   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	rhoReal1   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	
	fftw_Ctmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	fftw_Rtmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	gradient  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	gradientF  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	cplxTmp   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	quadTerm  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	cubTerm   = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	quarTerm  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	
	phi_F     = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	phi_R     = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	
	LaplaceR    = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	LaplaceF    = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	NonlinearR  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	NonlinearF  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofs);
	
	Gsquare   = (double *)malloc(sizeof(double) *cplxDofs);
	epmckt2   = (double *)malloc(sizeof(double) *cplxDofs);
    ikt2      = (double *)malloc(sizeof(double) *cplxDofs);


	for(int i = 0; i < cplxDofs; i++)
	{
		setCplxZero(rhoCplx[i]);
		setCplxZero(rhoReal[i]);
        setCplxZero(rhoReal1[i]);
        setCplxZero(rhoCplx[i]);
		setCplxZero(phi_F[i]);
		setCplxZero(gradient[i]);
		setCplxZero(gradientF[i]);
		setCplxZero(cplxTmp[i]);
		Gsquare[i] = 0.0;
		epmckt2[i] = 0.0;
        ikt2[i]    = 0.0;
		for(int j = 0; j < DimPhy; j++) projPlane[i][j] = 0.0;
		for(int k = 0; k < DimCpt; k++) indKspace[i][k] = 0;
	}

	Planc2cFord = fftw_plan_dft(DimCpt, NCpt, fftw_Rtmp, fftw_Ctmp, FFTW_FORWARD,  FFTW_MEASURE);  // real to cplx
	Planc2cBack = fftw_plan_dft(DimCpt, NCpt, fftw_Ctmp, fftw_Rtmp, FFTW_BACKWARD, FFTW_MEASURE);  // cplx to real 
}

void getProjPlane()
{
	double mnt;
	for(int i = 0; i < cplxDofs; i++)
	{
		for(int j = 0; j < DimPhy; j++)
		{
			mnt = 0.0;
            for(int k=0;k<DimCpt;k++)
            {
                mnt = mnt + rcpBox[j][k]*indKspace[i][k];
                //printf("indKspace[%d][%d] = %d ",i,k,indKspace[i][k]);
            }
            projPlane[i][j] = mnt;
            //printf("indKspace[%d][%d] = %.3e ",i,j,projPlane[i][j]);
		}
        //printf("\n");
	}
}

void getGsquare()
{
	double s;

	for(int i = 0; i < cplxDofs; i++)
	{
        s = 0;
        for(int j=0;j<DimPhy;j++)
        {
            s = s + pow(projPlane[i][j],2);
        }
        Gsquare[i] = s;
        //printf("Gsquare[%d] = %.3e ",i,Gsquare[i]);
        //printf("\n");
	}
}

void get_epmckt2()
{
    for(int i=0;i<cplxDofs;i++)
    {
        epmckt2[i] = (PARA_Q0*PARA_Q0-Gsquare[i])*(PARA_Q1*PARA_Q1-Gsquare[i]);
        epmckt2[i] = PARA_XI*pow(epmckt2[i],2);
        //printf("epmckt2[%d] = %.3e \n",i,epmckt2[i]);
    }
}

void get_ikt2()
{
    for(int i=0;i<cplxDofs;i++)
    {
        ikt2[i] = 1.0/( epmckt2[i] + abs(PARA_tau) + 0.02);
    }
}


void get_initial_value(fftw_complex *v,int choice, double angle)
{
    p = fftw_plan_dft(DimCpt,NCpt,v,fftw_Ctmp,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    FFT_dot_constant(fftw_Ctmp,cplxDofs,1.0/cplxDofs);

    printf("\nchoice = %d \t angle = %.3e ",choice,angle);
    //转化为弧度
    angle = (angle/180.0)*PI;
    int N = NCpt[1];
    int n = 0;
    int Kindex[8][24] = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23},
                         {0,2,4,6,8,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {12,14,16,18,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,1,3,5,6,7,9,11,12,17,18,23,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,1,6,7,12,14,15,16,18,20,21,22,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,1,6,7,12,18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
                        };

    //Kindex = (int *)malloc(sizeof(int) *24);
    if(choice==1)
    {
        cout << "\n--------------------QC----------------------------\n" << endl;
        n = 24;
    }
    else if(choice == 2)
    {
        cout << "\n--------------------C6----------------------------\n" << endl;
        n = 6;
    }
    else if(choice == 3)
    {
    	cout << "\n--------------------CS6----------------------------\n" << endl;
        n = 6;
    }
    else if(choice == 4)
    {
    	cout << "\n--------------------LQ----------------------------\n" << endl;
        n = 12;
    }
    else if(choice == 5)
    {
    	cout << "\n--------------------LSQ----------------------------\n" << endl;
        n = 12;
    }
    else if(choice == 6)
    {
        cout << "\n--------------------T6----------------------------\n" << endl;
        n = 6;
    }
    else if(choice == 7)
    {
    	cout << "\n--------------------HG----------------------------\n" << endl;
        n = 0;
    }
    else if(choice == 8)
    {
    	cout << "\n--------------------LAM----------------------------\n" << endl;
        n = 2;
    }
    
    for(int i=0;i<cplxDofs;i++)
    {
        fftw_Ctmp[i][0] = 0.0;
        fftw_Ctmp[i][1] = 0.0;
    }

    int ki;
    int k1,k2;
    double kx,ky;
    for(int i=0;i<n;i++)
    {
        ki = Kindex[choice-1][i];
        if(ki>11)
        {
            ki = ki-24;
            kx = PARA_Q1*cos(ki*2.0*PI/12.0 + PI/12.0 + angle);
            ky = PARA_Q1*sin(ki*2.0*PI/12.0 + PI/12.0 + angle);
        }
        else
        {
            kx = cos(ki*2.0*PI/12.0 + angle);
            ky = sin(ki*2.0*PI/12.0 + angle);
        }
        k1 = round(kx*L);
        k2 = round(ky*L);
        
        if(k1<0)
            k1 = N+k1;
        
        if(k2<0)
            k2 = N + k2;
        fftw_Ctmp[k1*N+k2][0] = 0.05;
    }

    q = fftw_plan_dft(DimCpt,NCpt,fftw_Ctmp,v,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(q);
}

















