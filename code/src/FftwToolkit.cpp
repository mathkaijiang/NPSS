#include "Head.h"
#include "Data.h"
#include "Initialization.h"
#include "FftwToolkit.h"
#include "BasicOperators.h"
#include "Mytimer.h"

void getIndex(int **kspace, int dim) 
{
	int *k = (int *)malloc(sizeof(int) *cplxDofs);
	for(int i = 0; i < cplxDofs; i++)
    {   
        if(i<NCpt[0]/2)
            k[i] = i;
        else
            k[i] = i-NCpt[0];
    }

	for(int k1=0;k1<NCpt[0];k1++)
    {
        for(int k2=0;k2<NCpt[0];k2++)
        {
            kspace[k1*NCpt[0]+k2][0] = k[k1];
            kspace[k1*NCpt[0]+k2][1] = k[k2];
        }
    }
	free(k);
}

int getIndex1D(int *k, int *Ndof, int dim)
{
	int index=0,hi;
    for(int i=0;i<dim;i++)
    {
        hi = 0;
        if(k[i]< 0)
            hi = k[i] + Ndof[i];
        else
            hi = k[i];

        index = index + hi*pow(Ndof[i],dim-1-i); 
    }
	return index;
}

void convolution(fftw_complex *src, fftw_complex *orig, int order)
{
	memcpy(fftw_Ctmp, src, sizeof(fftw_complex)*cplxDofs);
	fftw_execute(Planc2cBack);   // cplx to real
	for(int i = 0; i < cplxDofs; i++)
	{
		fftw_Rtmp[i][0] = std::pow(fftw_Rtmp[i][0], order);
		fftw_Rtmp[i][1] = std::pow(fftw_Rtmp[i][1], order);
	}
	fftw_execute(Planc2cFord);  // real to cplx
	FuncsLinear1Cplx(fftw_Ctmp, cplxDofs, 1.0/cplxDofs, fftw_Ctmp);
	memcpy(orig, fftw_Ctmp, sizeof(fftw_complex)*cplxDofs);
}


void pad(fftw_complex *in, fftw_complex *out, int *Ndof, int *ExpandNdof, int dim)
{
	int k, Expandk;
	for(int i = 0; i < cplxDofs; i++)
	{
		k = getIndex1D(indKspace[i], Ndof, dim);
		Expandk = getIndex1D(indKspace[i], ExpandNdof, dim);
		memcpy(out[Expandk], in[k], sizeof(fftw_complex)*1);
//        out[Expandk][0] = in[k][0]; out[Expandk][1] = in[k][1];
	}
}
