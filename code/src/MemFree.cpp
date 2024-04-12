#include "Data.h"
#include "MemFree.h"

void memReleaser()
{
	fftw_free(rhoCplx);
	
	fftw_free(rhoReal);
    fftw_free(rhoCplx1);
    fftw_free(rhoReal1);

	fftw_free(fftw_Ctmp);
	fftw_free(fftw_Rtmp);
	fftw_free(gradient);
	fftw_free(gradientF);
	fftw_free(cplxTmp);
	fftw_free(quadTerm);
	fftw_free(cubTerm);
	fftw_free(quarTerm);
	fftw_free(phi_F);
	fftw_free(phi_R);
	
	free(Gsquare);
	free(epmckt2);
    free(ikt2);
	
	fftw_free(LaplaceR);
	fftw_free(LaplaceF);
	fftw_free(NonlinearR);
	fftw_free(NonlinearF);
	
	for(int i = 0; i < DimCpt; i++) 
	{
		free(dirBox[i]); 
		free(rcpBox[i]);
	}
	free(dirBox);
	free(rcpBox);
	
	for(int i = 0; i < DimPhy; i++) free(ProjMatrix[i]);
	free(ProjMatrix);
	
	//printf("\n go go go \n ");
	for(int i = 0; i < cplxDofs; i++)
	{
		free(indKspace[i]);
		free(projPlane[i]);
		//printf("\n go go go \n ");
	}
	
	//printf("\n go go go \n ");
	
	free(indKspace);
	free(projPlane);
	//printf("\n go go go \n ");
	
	fftw_destroy_plan(Planc2cFord);
	fftw_destroy_plan(Planc2cBack);
	fftw_destroy_plan(p);
	fftw_destroy_plan(q);
	//printf("\n go go go \n ");
}
