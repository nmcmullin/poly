#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cmath>
#include<vector>
#include<string>
#include <algorithm>
#include <numeric>
#include <sstream>
#include "CFD.h"
#include "NavierStokesCalc.h"
#include "FreeSurfaceCalc.h"
#include "CFDRunClass.h"
#include "flag.h"
#define PI 3.14159
#define g 9.81
using namespace std;
int main()
{
	ofstream foutu("InitialU.dat",ios::binary);
	ofstream foutv("InitialV.dat",ios::binary);
	ofstream foutp("InitialP.dat",ios::binary);
	ofstream foutm("mask.dat",ios::binary);
	ofstream foutx("InitialX.dat",ios::binary);
	ofstream fouty("InitialY.dat",ios::binary);
	ofstream foutc("InitialC.dat",ios::binary);

int ncolsT=82;
int nrowsT=82;
float dx=0.5;
float dy=0.25;
float x1;
float y1;
float C1;
vector<int> z(nrowsT*ncolsT);
vector<int> x;
vector<int> y;
vector<float> C(nrowsT*ncolsT);
vector<int> maskT((nrowsT+2)*(ncolsT+2));
/*for(int j=1;j<=(3*(nrowsT-2)/4);j++)
{
		for(int k=1;k<=3;k++)
		{
			for(int i=1;i<=((ncolsT-2)/2);i++)
			{			
				for(int l=1;l<=3;l++)
				{
				x1=((i-1)*dx)+(l*(dx/4));
				//x.push_back(x1);
				foutx<<x1<<" ";
				y1=((j-1)*dy)+(k*(dy/4));
				//y.push_back(y1);
				fouty<<y1<<" ";
				}
			}
			foutx<<endl;
			fouty<<endl;
		}
}
for(int j=1;j<=(3*(nrowsT-2)/4);j++)
{
	for(int i=1;i<=((ncolsT-2)/2);i++)
	{
	C[j*ncolsT+i]=1.0;
	}
	C[j*ncolsT+((ncolsT-2)/2)+1]=0.5;
}
for(int i=1;i<=((ncolsT-2)/2)+1;i++)
{
	C[((3*(nrowsT-2)/4)+1)*ncolsT+i]=0.5;
}*/
for(int j=1;j<=(nrowsT-2);j++)
{
	for(int i=1;i<=(ncolsT-2);i++)
	{
		C[j*ncolsT+i]=0.0;
	}

}
for(int j=2;j<=(nrowsT-1);j++)
	{
		for(int i=2;i<=(ncolsT-1);i++)
		{
			maskT[j*ncolsT+i]=0;
		}
	}
for(int j=0;j<=(nrowsT+1);j++)
	{
		maskT[j*ncolsT+0]=1;
		maskT[j*ncolsT+1]=1;
		maskT[j*ncolsT+ncolsT]=1;
		maskT[j*ncolsT+ncolsT+1]=1;
	}
for(int i=2;i<=(ncolsT-1);i++)
	{
		maskT[0*ncolsT+i]=1;
		maskT[1*ncolsT+i]=1;
		maskT[nrowsT*ncolsT+i]=1;
		maskT[(nrowsT+1)*ncolsT+i]=1;
	}

for(int j=0;j<=(nrowsT-1);j++)
	{
		for(int i=0;i<=(ncolsT-1);i++)
		{	
			z[j*ncolsT+i]=0;
			foutu<<z[j*ncolsT+i]<<" ";
			foutv<<z[j*ncolsT+i]<<" ";
			foutp<<z[j*ncolsT+i]<<" ";
			foutc<<C[j*ncolsT+i]<<" ";

		}
		foutu<<endl;
		foutv<<endl;
		foutp<<endl;
		foutc<<endl;
	}

for(int j=0;j<=(nrowsT+1);j++)
	{
		for(int i=0;i<=(ncolsT+1);i++)
		{	
			foutm<<maskT[j*ncolsT+i]<<" ";
		}
		foutm<<endl;
	}
	float lambda=20;
	float rho=1000;
	float dynvis=1.5e-3;
	float omega=10;
	float k=2*PI/lambda;
	float velocityscale=omega/k;
	float lengthscale=2*2*lambda;
	float Re=(rho/dynvis)*velocityscale*lengthscale;
																														//if SMAC?=false, npart=0
					//ncols, nrows, dx, dy Re, tmax, maxNumIt, veleps, peps, gx, gy, w, tsf, gsf, GausSeidel?, Conjugate Gradient? SMAC? npart, SLIC?, FLAP?, YVOF?,voftol
	CFDRunClass CFDRC1(82,82,0.5,0.25,Re,rho,dynvis,3*PI/20,500,1.0e-7,1.0e-7,0.0,-g,1.7,0.5,1.5,true,false,false,0,false,false,true,0.0);//npart27648 
	CFDRC1.Run();
	return 0;
}
