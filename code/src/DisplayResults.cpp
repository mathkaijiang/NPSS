#include "Data.h"
#include "Initialization.h"
#include "FftwToolkit.h"
#include "BasicOperators.h"
#include "Mytimer.h"
#include "DisplayResults.h"

void dispDensity(fftw_complex *src, int step)
{
	mytimer_t timer;
	timer.reset();
	timer.start();

	FILE *densityfile;
	char densityName[100];
	sprintf(densityName, "./result/density_%s_%d.dat",model_initial_value,realDofs);
	densityfile = fopen(densityName, "w");

	if(DimCpt == 3 && DimPhy == 3)
	{
		int *ExpandNCpt = (int *)malloc(sizeof(int)*DimCpt);
		for(int i = 0; i < DimCpt; i++) ExpandNCpt[i] = 48;
		int Expandof = 1;
		for(int i = 0; i < DimCpt; i++) Expandof *= ExpandNCpt[i];

		fftw_complex *ExpandRho = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*Expandof);
		fftw_complex *padRho = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*Expandof);
		for(int i = 0; i < Expandof; i++) setCplxZero(ExpandRho[i]);
		fftw_plan Pc2r = fftw_plan_dft(DimCpt, ExpandNCpt, ExpandRho, padRho, FFTW_BACKWARD, FFTW_MEASURE);  // cplx to real 

		pad(src, ExpandRho, NCpt, ExpandNCpt, DimCpt);
		fftw_execute(Pc2r);

		int n0 = ExpandNCpt[0];
		int n1 = ExpandNCpt[1];
		int n2 = ExpandNCpt[2];

		for(int i=0;i<n0;i++)	fprintf(densityfile,"%.10e\t",i*dirBox[0][0]/(double)n0);
		fprintf(densityfile,"\n");

		for(int i=0;i<n1;i++)	fprintf(densityfile,"%.10e\t",i*dirBox[1][1]/(double)n1);
		fprintf(densityfile,"\n");

		for(int i=0;i<n2;i++)	fprintf(densityfile,"%.10e\t",i*dirBox[2][2]/(double)n2);
		fprintf(densityfile,"\n");

		int M=0;	
		for(int i=0;i<n0;i++)
		{
			for(int j=0;j<n1;j++)
			{
				for(int k=0;k<n2;k++)
				{
					fprintf(densityfile,"%.10e ", padRho[M][0]);
					M=M+1;
				}
				fprintf(densityfile,"\n");
			}
			fprintf(densityfile,"\n");
		}
		free(ExpandNCpt);
		fftw_free(ExpandRho);
		fftw_free(padRho);
		fftw_destroy_plan(Pc2r);
		fclose(densityfile);
	}

	if(DimCpt == 4 && DimPhy == 2)
	{
		double dy = 0.7;
		double dz = 0.7;
		double enlarge = 30;
		double sizey = enlarge*dirBox[0][0];
		double sizez = enlarge*dirBox[0][0];
		int ny, nz;  // 实空间的离散点
		ny = ceil(sizey/dy);
		nz = ceil(sizez/dz);
		double tol = 1.0e-5;
		
		for(int ky = 0; ky < ny+1; ky++)
		{
			for(int kz = 0; kz < nz+1; kz++)
			{
				double rho = 0.0;
				double tmpphase;

				for(int i = 0; i < cplxDofs; i ++)
				{
					double elm = std::sqrt(src[i][0]*src[i][0]+src[i][1]*src[i][1]);
					if(elm > tol)
					{
						tmpphase = projPlane[i][0]*dy*ky+projPlane[i][1]*dz*kz;
						rho += (src[i][0]*cos(tmpphase)-src[i][1]*sin(tmpphase));
					}
				}
				fprintf(densityfile, "%.15e\t", rho);
			}
			fprintf(densityfile, "\n");
		}
	}
	timer.pause();
	printf("\n===> Output plotted data, Step %d: %f seconds\n\n", step, timer.get_current_time());
}

bool myComp(mySortVec a, mySortVec b)
{
	return (a.Data[0]*a.Data[0]+a.Data[1]*a.Data[1] > b.Data[0]*b.Data[0]+b.Data[1]*b.Data[1]);
}

bool scaleComp(mySortWave a, mySortWave b)
{
	return (a.Data > b.Data);
}

void dispCoeffFourier(fftw_complex *src, int step)
{
	FILE *fFourCoeff;
	char fname[100];
	sprintf(fname, "./result/field.%d.dat", step);
	fFourCoeff = fopen(fname, "w");

	mySortVec *myVector = (mySortVec*)malloc(sizeof(mySortVec)*cplxDofs);
	for(int k = 0; k < cplxDofs; k++)
	{
		myVector[k].Data  = (double *)malloc(sizeof(double)*2);
		myVector[k].Index = (int *)malloc(sizeof(int)*DimCpt);
		for(int i = 0; i < DimCpt; i++)
			myVector[k].Index[i] = indKspace[k][i];
		for(int j = 0; j < 2; j++)
			myVector[k].Data[j] = src[k][j];
	}

	std::sort(myVector, myVector+cplxDofs, myComp);
	for(int k = 0; k < cplxDofs; k++)
	{
		for(int i = 0; i < DimCpt; i++)
		{
			fprintf(fFourCoeff, "%d\t", myVector[k].Index[i]);
		}
		fprintf(fFourCoeff, "%e\t%e\n", myVector[k].Data[0], myVector[k].Data[1]);
	}

	for(int i = 0; i < cplxDofs; i++)
	{
		free(myVector[i].Index);
		free(myVector[i].Data);
	}
	free(myVector);
	fclose(fFourCoeff);
}

void dispPlaneWave(fftw_complex *src)
{
	FILE *fprojPlane;
	fprojPlane = fopen("./result/dispPlaneWave.dat", "w");

	mySortWave *myPlaneWave = (mySortWave*)malloc(sizeof(mySortWave)*cplxDofs);
	for(int k = 0; k < cplxDofs; k++)
	{
		myPlaneWave[k].Wave = (double *)malloc(sizeof(double)*DimPhy);
		for(int i = 0; i < DimPhy; i++)
			myPlaneWave[k].Wave[i] = projPlane[k][i];
		myPlaneWave[k].Data = std::sqrt(src[k][0]*src[k][0]+src[k][1]*src[k][1]);
	}

	std::sort(myPlaneWave, myPlaneWave+cplxDofs, scaleComp);

	for(int i = 0; i < cplxDofs; i ++)
	{
		for(int j = 0; j < DimPhy; j++)
			fprintf(fprojPlane, "%f\t ", myPlaneWave[i].Wave[j]);
		fprintf(fprojPlane, "%e\n ", myPlaneWave[i].Data);
	}
	
	for(int i = 0; i < cplxDofs; i++)
	{
		free(myPlaneWave[i].Wave);
	}
	free(myPlaneWave);
	fclose(fprojPlane);
}
