#include "Head.h"
#include "DualBox.h"

void getRecipLattice(double **dBox, double **rBox, int dim)
{
	if(dim==2)
	{
		double det = dBox[0][0]*dBox[1][1]-dBox[0][1]*dBox[1][0];
		rBox[0][0] =  2*PI*dBox[1][1] / det;
		rBox[1][0] = -2*PI*dBox[0][1] / det;
		rBox[0][1] = -2*PI*dBox[1][0] / det;
		rBox[1][1] =  2*PI*dBox[0][0] / det;
	}
	if(dim == 3)
	{
		double aone[3], atwo[3], athree[3];
		double bone[3], btwo[3], bthree[3];
		double cone[3], ctwo[3], cthree[3];
		double volume[3];

		for(int i = 0; i < 3; i++)
		{
			aone[i]   = dBox[0][i];
			atwo[i]   = dBox[1][i];
			athree[i] = dBox[2][i];
		}

		cone[0]=atwo[1]*athree[2]-atwo[2]*athree[1];
		cone[1]=atwo[2]*athree[0]-atwo[0]*athree[2];
		cone[2]=atwo[0]*athree[1]-atwo[1]*athree[0];
		volume[0]=aone[0]*cone[0]+aone[1]*cone[1]+aone[2]*cone[2];
		bone[0]=2*PI/volume[0]*cone[0];
		bone[1]=2*PI/volume[0]*cone[1];
		bone[2]=2*PI/volume[0]*cone[2];    

		ctwo[0]=athree[1]*aone[2]-athree[2]*aone[1];
		ctwo[1]=athree[2]*aone[0]-athree[0]*aone[2];
		ctwo[2]=athree[0]*aone[1]-athree[1]*aone[0];
		volume[1]=atwo[0]*ctwo[0]+atwo[1]*ctwo[1]+atwo[2]*ctwo[2];
		btwo[0]=2*PI/volume[1]*ctwo[0];
		btwo[1]=2*PI/volume[1]*ctwo[1];
		btwo[2]=2*PI/volume[1]*ctwo[2];  
			
		cthree[0]=aone[1]*atwo[2]-aone[2]*atwo[1];
		cthree[1]=aone[2]*atwo[0]-aone[0]*atwo[2];
		cthree[2]=aone[0]*atwo[1]-aone[1]*atwo[0];
		volume[2]=athree[0]*cthree[0]+athree[1]*cthree[1]+athree[2]*cthree[2];
		bthree[0]=2*PI/volume[2]*cthree[0];
		bthree[1]=2*PI/volume[2]*cthree[1];
		bthree[2]=2*PI/volume[2]*cthree[2];   
		for(int i = 0; i < 3; i++)
		{
			rBox[0][i] = bone[i];
			rBox[1][i] = btwo[i]; 
			rBox[2][i] = bthree[i];
		}
	}
	if(dim == 4)
	{
		rBox[0][0] = 2.0*PI / dBox[0][0];
		rBox[1][1] = 2.0*PI / dBox[1][1];
		rBox[2][2] = 2.0*PI / dBox[2][2];
		rBox[3][3] = 2.0*PI / dBox[3][3];
//	 超六面体
//        double det =  dBox[0][0]*dBox[1][1]*dBox[2][2]*dBox[3][3] - dBox[0][0]*dBox[1][1]*dBox[2][3]*dBox[3][2] 
//                    - dBox[0][0]*dBox[1][2]*dBox[2][1]*dBox[3][3] + dBox[0][0]*dBox[1][2]*dBox[2][3]*dBox[3][1] 
//                    + dBox[0][0]*dBox[1][3]*dBox[2][1]*dBox[3][2] - dBox[0][0]*dBox[1][3]*dBox[2][2]*dBox[3][1] 
//                    - dBox[0][1]*dBox[1][0]*dBox[2][2]*dBox[3][3] + dBox[0][1]*dBox[1][0]*dBox[2][3]*dBox[3][2] 
//                    + dBox[0][1]*dBox[1][2]*dBox[2][0]*dBox[3][3] - dBox[0][1]*dBox[1][2]*dBox[2][3]*dBox[3][0] 
//                    - dBox[0][1]*dBox[1][3]*dBox[2][0]*dBox[3][2] + dBox[0][1]*dBox[1][3]*dBox[2][2]*dBox[3][0] 
//                    + dBox[0][2]*dBox[1][0]*dBox[2][1]*dBox[3][3] - dBox[0][2]*dBox[1][0]*dBox[2][3]*dBox[3][1] 
//                    - dBox[0][2]*dBox[1][1]*dBox[2][0]*dBox[3][3] + dBox[0][2]*dBox[1][1]*dBox[2][3]*dBox[3][0] 
//                    + dBox[0][2]*dBox[1][3]*dBox[2][0]*dBox[3][1] - dBox[0][2]*dBox[1][3]*dBox[2][1]*dBox[3][0] 
//                    - dBox[0][3]*dBox[1][0]*dBox[2][1]*dBox[3][2] + dBox[0][3]*dBox[1][0]*dBox[2][2]*dBox[3][1] 
//                    + dBox[0][3]*dBox[1][1]*dBox[2][0]*dBox[3][2] - dBox[0][3]*dBox[1][1]*dBox[2][2]*dBox[3][0] 
//                    - dBox[0][3]*dBox[1][2]*dBox[2][0]*dBox[3][1] + dBox[0][3]*dBox[1][2]*dBox[2][1]*dBox[3][0];
//        rBox[0][0] = 2*PI* (dBox[1][1]*dBox[2][2]*dBox[3][3] - dBox[1][1]*dBox[2][3]*dBox[3][2] 
//                         - dBox[1][2]*dBox[2][1]*dBox[3][3] + dBox[1][2]*dBox[2][3]*dBox[3][1] 
//                         + dBox[1][3]*dBox[2][1]*dBox[3][2] - dBox[1][3]*dBox[2][2]*dBox[3][1]) / det;
//        rBox[0][1] = -(2*PI*(dBox[1][0]*dBox[2][2]*dBox[3][3] - dBox[1][0]*dBox[2][3]*dBox[3][2] 
//                        - dBox[1][2]*dBox[2][0]*dBox[3][3] + dBox[1][2]*dBox[2][3]*dBox[3][0] 
//                        + dBox[1][3]*dBox[2][0]*dBox[3][2] - dBox[1][3]*dBox[2][2]*dBox[3][0])) / det;
//        rBox[0][2] = 2*PI*(dBox[1][0]*dBox[2][1]*dBox[3][3] - dBox[1][0]*dBox[3][1]*dBox[2][3] 
//                        - dBox[1][1]*dBox[2][0]*dBox[3][3] + dBox[1][1]*dBox[3][0]*dBox[2][3] 
//                        + dBox[2][0]*dBox[1][3]*dBox[3][1] - dBox[2][1]*dBox[3][0]*dBox[1][3]) / det;	
//        rBox[0][3] = -2*PI*(dBox[1][0]*dBox[2][1]*dBox[3][2] - dBox[1][0]*dBox[2][2]*dBox[3][1] 
//                        - dBox[1][1]*dBox[2][0]*dBox[3][2] + dBox[1][1]*dBox[3][0]*dBox[2][2] 
//                        + dBox[2][0]*dBox[1][2]*dBox[3][1] - dBox[1][2]*dBox[2][1]*dBox[3][0]) / det;
//        rBox[1][0] = -(2*PI*(dBox[0][1]*dBox[2][2]*dBox[3][3] - dBox[0][1]*dBox[2][3]*dBox[3][2] 
//                        - dBox[0][2]*dBox[2][1]*dBox[3][3] + dBox[0][2]*dBox[3][1]*dBox[2][3] 
//                        + dBox[0][3]*dBox[2][1]*dBox[3][2] - dBox[0][3]*dBox[2][2]*dBox[3][1])) / det;
//        rBox[1][1] = (2*PI*(dBox[0][0]*dBox[2][2]*dBox[3][3] - dBox[0][0]*dBox[2][3]*dBox[3][2] 
//                        - dBox[0][2]*dBox[2][0]*dBox[3][3] + dBox[0][2]*dBox[3][0]*dBox[2][3] 
//                        + dBox[2][0]*dBox[0][3]*dBox[3][2] - dBox[0][3]*dBox[3][0]*dBox[2][2])) / det;
//        rBox[1][2] = -(2*PI*(dBox[0][0]*dBox[2][1]*dBox[3][3] - dBox[0][0]*dBox[3][1]*dBox[2][3] 
//                        - dBox[0][1]*dBox[2][0]*dBox[3][3] + dBox[0][1]*dBox[3][0]*dBox[2][3] 
//                        + dBox[2][0]*dBox[0][3]*dBox[3][1] - dBox[0][3]*dBox[2][1]*dBox[3][0])) / det;
//        rBox[1][3] = (2*PI*(dBox[0][0]*dBox[2][1]*dBox[3][2] - dBox[0][0]*dBox[2][2]*dBox[3][1] 
//                        - dBox[0][1]*dBox[2][0]*dBox[3][2] + dBox[0][1]*dBox[3][0]*dBox[2][2] 
//                        + dBox[0][2]*dBox[2][0]*dBox[3][1] - dBox[0][2]*dBox[2][1]*dBox[3][0])) / det;
//        rBox[2][0] = (2*PI*(dBox[0][1]*dBox[1][2]*dBox[3][3] - dBox[0][1]*dBox[1][3]*dBox[3][2]
//                        - dBox[0][2]*dBox[1][1]*dBox[3][3] + dBox[0][2]*dBox[1][3]*dBox[3][1] 
//                        + dBox[1][1]*dBox[0][3]*dBox[3][2] - dBox[0][3]*dBox[1][2]*dBox[3][1])) / det;
//        rBox[2][1] = -(2*PI*(dBox[0][0]*dBox[1][2]*dBox[3][3] - dBox[0][0]*dBox[1][3]*dBox[3][2] 
//                        - dBox[1][0]*dBox[0][2]*dBox[3][3] + dBox[1][0]*dBox[0][3]*dBox[3][2] 
//                        + dBox[0][2]*dBox[3][0]*dBox[1][3] - dBox[0][3]*dBox[1][2]*dBox[3][0])) / det;
//        rBox[2][2] = (2*PI*(dBox[0][0]*dBox[1][1]*dBox[3][3] - dBox[0][0]*dBox[1][3]*dBox[3][1] 
//                        - dBox[0][1]*dBox[1][0]*dBox[3][3] + dBox[0][1]*dBox[3][0]*dBox[1][3] 
//                        + dBox[1][0]*dBox[0][3]*dBox[3][1] - dBox[1][1]*dBox[0][3]*dBox[3][0])) / det;
//        rBox[2][3] = -(2*PI*(dBox[0][0]*dBox[1][1]*dBox[3][2] - dBox[0][0]*dBox[1][2]*dBox[3][1] 
//                        - dBox[0][1]*dBox[1][0]*dBox[3][2] + dBox[0][1]*dBox[1][2]*dBox[3][0] 
//                        + dBox[1][0]*dBox[0][2]*dBox[3][1] - dBox[0][2]*dBox[1][1]*dBox[3][0])) / det;
//        rBox[3][0] = -(2*PI*(dBox[0][1]*dBox[1][2]*dBox[2][3] - dBox[0][1]*dBox[1][3]*dBox[2][2] 
//                        - dBox[0][2]*dBox[1][1]*dBox[2][3] + dBox[0][2]*dBox[2][1]*dBox[1][3] 
//                        + dBox[1][1]*dBox[0][3]*dBox[2][2] - dBox[0][3]*dBox[1][2]*dBox[2][1])) / det;
//        rBox[3][1] = (2*PI*(dBox[0][0]*dBox[1][2]*dBox[2][3] - dBox[0][0]*dBox[1][3]*dBox[2][2] 
//                        - dBox[1][0]*dBox[0][2]*dBox[2][3] + dBox[1][0]*dBox[0][3]*dBox[2][2] 
//                        + dBox[0][2]*dBox[2][0]*dBox[1][3] - dBox[2][0]*dBox[0][3]*dBox[1][2])) / det;
//        rBox[3][2] = -(2*PI*(dBox[0][0]*dBox[1][1]*dBox[2][3] - dBox[0][0]*dBox[2][1]*dBox[1][3] 
//                        - dBox[0][1]*dBox[1][0]*dBox[2][3] + dBox[0][1]*dBox[2][0]*dBox[1][3] 
//                        + dBox[1][0]*dBox[0][3]*dBox[2][1] - dBox[1][1]*dBox[2][0]*dBox[0][3])) / det;
//        rBox[3][3] = (2*PI*(dBox[0][0]*dBox[1][1]*dBox[2][2] - dBox[0][0]*dBox[1][2]*dBox[2][1] 
//                        - dBox[0][1]*dBox[1][0]*dBox[2][2] + dBox[0][1]*dBox[2][0]*dBox[1][2] 
//                        + dBox[1][0]*dBox[0][2]*dBox[2][1] - dBox[0][2]*dBox[1][1]*dBox[2][0])) / det;
	}
}
