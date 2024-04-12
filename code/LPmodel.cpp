#include "Data.h"
#include "Head.h"
#include "Initialization.h"
#include "FftwToolkit.h"
#include "BasicOperators.h"
#include "DisplayResults.h"
#include "LPsrc.h"
#include "NPSS.h"
#include "HISD.h"
#include "Mytimer.h"
#include "MemFree.h"
#include<ctime>

void Save_phi(fftw_complex *x, char *filename);
void Save_Option_V(Option *option, char *filename);
void load_Option_V(Option *option, char *filename);
void load_phi(fftw_complex *x, char *filename);
using namespace std;

int main()
{
	printf("\n\n ====================================   START PROGRAM  ======================================\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	double hamilton = 0.0,ene=0,res = 0,res1 = 0;
	double TOL = 1.0e-20;
	char filename[100];
	clock_t startTime,endTime;

	initialize_LP();
    //hamilton = gradientflow(rhoReal, TOL,2000);
	sprintf(filename, "./result/%s","x0.txt");	
	load_phi(rhoReal,filename);

    ene = ene_cam(rhoReal);
	//hamilton = ene;
    ngrad_cam(rhoReal,gradient);
    res = nrmp(gradient);
    printf("\n ene: %.3e \n",ene);
    printf("\n res: %.13e \n ",res);

   
    
    option1 = (Option *)malloc(sizeof(Option));
    hiosd_lobpcg_initialization(option1,5);
    printf("\n hhh \n");

	sprintf(filename, "./result/%s","V0.txt");
	load_Option_V(option1,  filename);	
	//for(int i=0;i<20;i++)
	//{
	//	hiosd_lobpcg_inip(rhoReal,option1);
    //}

	/* ---------------------------save--------------------------------------*/
	//sprintf(filename, "./result/%s","x0.txt");
	//Save_phi(rhoReal,filename);
	//sprintf(filename, "./result/%s","V0.txt");
	//Save_Option_V(option1, filename);	
	/*----------------------------------------------------------------------*/

	/*---------------------------NPSS----------------------------------------*/
    startTime = clock();
	option2 = (Option *)malloc(sizeof(Option));
    hiosd_lobpcg_initialization(option2,1);
    critical_point();
    
	//sprintf(filename, "./result/%s","x1_GOSD.txt");	
	//load_phi(rhoReal1,filename);

    endTime = clock();

	for(int i = 0;i<2;i++)
	{
		hiosd_lobpcg_inip(rhoReal1,option1);
	}

    printf("\n--------------------end--------------------\n");
	/*----------------------------------------------------------------------*/

	/* ---------------------------save--------------------------------------*/
	sprintf(filename, "./result/%s","x1.txt");
	Save_phi(rhoReal1,filename);
	sprintf(filename, "./result/%s","V1_NPSS.txt");
	Save_Option_V(option1, filename);	
	/*----------------------------------------------------------------------*/
	
	ene = ene_cam(rhoReal1);
	//hamilton = ene;
    ngrad_cam(rhoReal1,gradient);
    res = nrmp(gradient);
    printf("\n ene: %.3e \n",ene);
    printf("\n res: %.13e \n ",res);


    for(int i=0;i<option1->k;i++)
    {
        fftw_free(option1->V[i]);
    }
    
    for(int i=0;i<option2->k;i++)
    {
        fftw_free(option2->V[i]);
    }
    
    fftw_free(option1->V);
    fftw_free(option2->V);
    free(option1);
    free(option2);

    printf("\n ene: %.3e \n",ene);
    printf("\n res: %.13e \n ",res);

	memReleaser();
	timer.pause();
	
	double csd_time = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	printf("计算鞍点 CPU 占用的时间：%f seconds\n",csd_time);
	printf("\n\n\n\t\t time cost of program : %f seconds, hamilton = %.16e\n", timer.get_current_time(), hamilton);
	printf("\n\n ======================================   END PROGRAM  ======================================\n\n");
	return 1;
}

void Save_phi(fftw_complex *x, char *filename)
{
	FILE *fp;
    fp = fopen(filename,"w");
    for(int i=0;i<cplxDofs;i++)
    {
        fprintf(fp,"%.16lf\n",x[i][0]);
    }
    fclose(fp);
}

void Save_Option_V(Option *option, char *filename)
{
	FILE *fp;
	fp = fopen(filename,"w");
    for(int i=0;i<option->k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            fprintf(fp,"%.16lf\n",option->V[i][j][0]);
        }
        //fprintf(fp,"\n");
    }
    fclose(fp);
}

void load_phi(fftw_complex *x, char *filename)
{
	printf("\n -----------loading------------------");
	FILE *fp;
    fp = fopen(filename,"r");
    for(int i=0;i<cplxDofs;i++)
    {
        fscanf(fp,"%lf",&(x[i][0]));
    }
    fclose(fp);
}

void load_Option_V(Option *option, char *filename)
{
	printf("\n -----------loading------------------");
	FILE *fp;
	fp = fopen(filename,"r");
    for(int i=0;i<option->k;i++)
    {
        for(int j=0;j<cplxDofs;j++)
        {
            fscanf(fp,"%lf",&(option->V[i][j][0]));
        }
    }
    fclose(fp);
}
