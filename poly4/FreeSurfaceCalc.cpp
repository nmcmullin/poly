#include "FreeSurfaceCalc.h"
#include "CFD.h"
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
#include "flag.h"
#define PI 3.14159
//All ny calcs are using dx, might need to change angle finding algorithm/ make angle finding algorithm into separate function
//Constructor
FreeSurfaceCalc::FreeSurfaceCalc(int inncols,int innrows, float indx, float indy, float inRe, float inrho, float indynvis, bool inSMACbool,int innpart,bool inSLICbool,bool inFLAPbool,bool inYVOFbool, float invoftol):CFD(inncols,innrows,indx,indy,inRe,inrho,indynvis,inSMACbool,inSLICbool,inFLAPbool,inYVOFbool),npart(innpart),voftol(invoftol),initialdirectionbool(0)
{
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	Fb.push_back(0.0);
	Ft.push_back(0.0);
	Fl.push_back(0.0);
	Fr.push_back(0.0);
	}
}
}
FreeSurfaceCalc::~FreeSurfaceCalc()
{
}
//Methods
void FreeSurfaceCalc::RunFreeSurfaceCalc()
{
PeriodicVolumeFractionBC();
if(SMACbool==true)
{
//RunSMAC
UpdateParticleLocations();
UpdateFluidFlag();
SetEmptyCellstoZero();
}
if(SLICbool==true)
{
RunSLIC();
UpdateFluidFlag();
SetEmptyCellstoZero();
}
else if(FLAPbool==true)
{
RunFLAP();
UpdateFluidFlag();
SetEmptyCellstoZero();
}
else if(YVOFbool==true)
{
RunYVOF();
UpdateFluidFlag();
SetEmptyCellstoZero();
}
}
void FreeSurfaceCalc::UpdateParticleLocations()
{
int i,j;
	for(int n=0;n<npart;n++)
	{
		i=ceil(x[n]/dx);
		j=ceil(y[n]/dy);

			if(fl[j*ncols+i].F==1)
			{
			x[n]+=((u[j*ncols+i]+u[j*ncols+i-1])/2.0)*dt;
			y[n]+=((v[j*ncols+i]+v[(j-1)*ncols+i])/2.0)*dt;
			}
	}
}
void FreeSurfaceCalc::UpdateFluidFlag()
{
if(SMACbool==true)
{
for(int j=0;j<nrows;j++)
	{
	for(int i=0;i<ncols;i++)
		{	if(fl[j*ncols+i].C==true)
			{
			int n=0;
			bool flagmodified=false;
			
			while(n<npart && flagmodified==false)
			{
			int I=ceil(x[n]/dx);
			int J=ceil(y[n]/dy);
				if(i==I && j==J)
				{
				flagmodified=SetFluidFlagTrue(i,j);
				}
				n+=1;
			}
		if(flagmodified==false)
		{
			flagmodified=SetFluidFlagFalse(i,j);
		}
		}	
	}
	}
}
else if(SLICbool==true || FLAPbool==true || YVOFbool==true)
{
for(int j=0;j<nrows;j++)
	{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true)
	{
	if(C[j*ncols+i]>voftol)
	{
		fl[j*ncols+i].F=true;
		if(C[j*ncols+i]<(1.0-voftol))
		{
			fl[j*ncols+i].S=true;
		}
		else
		{
			fl[j*ncols+i].S=false;
		}
	}
	else
	{
		fl[j*ncols+i].F=false;
		fl[j*ncols+i].S=false;
	}
}
}
}
}
}
bool FreeSurfaceCalc::SetFluidFlagTrue(int &passI, int &passJ)
{
	fl[passJ*ncols+passI].F=true;
	return true;
}
bool FreeSurfaceCalc::SetFluidFlagFalse(int &passI, int &passJ)
{
	fl[passJ*ncols+passI].F=false;
	return true;
}
void FreeSurfaceCalc::SetEmptyCellstoZero()
 {
for(int j=0;j<nrows;j++)
{
		for(int i=0;i<ncols;i++)
		{
			if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==false)
			{
			u[j*ncols+i]=0.0;
			//u[j*ncols+i-1]=0.0;
			v[j*ncols+i]=0.0;
			//v[(j-1)*ncols+i]=0.0;
			p[j*ncols+i]=0.0;
			}
		}
}
}
void FreeSurfaceCalc::PrintFreeSurfaceData(int &myMval)
{
char bufferx[50];
int filenamex=sprintf(bufferx, "Xt000%d.dat", (myMval));
char buffery[50];
int filenamey=sprintf(buffery, "Yt000%d.dat", (myMval));
ofstream foutx(bufferx,ios::binary);
ofstream fouty(buffery,ios::binary);

for(int n=0;n<npart;n++)
{
		foutx<<x[n]<<" ";
		fouty<<y[n]<<" ";
}
}
void FreeSurfaceCalc::PeriodicVolumeFractionBC()
{
	for(int j=0;j<nrows;j++)
	{
	C[j*ncols+0]=C[j*ncols+ncols-2];
	C[j*ncols+ncols-1]=C[j*ncols+1];
	Fr[j*ncols+0]=Fr[j*ncols+ncols-2];
	Fl[j*ncols+ncols-2]=Fl[j*ncols+1];
	}
	for(int i=0;i<ncols;i++)
	{
	C[0*ncols+i]=C[1*ncols+i];
	C[(nrows-1)*ncols+i]=C[(nrows-2)*ncols+i];
	}
}
void FreeSurfaceCalc::VolumeFractionBC()
{
	for(int j=0;j<nrows;j++)
	{
	C[j*ncols+0]=C[j*ncols+ncols-2];
	C[j*ncols+ncols-1]=C[j*ncols+1];
	}
	for(int i=0;i<ncols;i++)
	{
	C[0*ncols+i]=C[1*ncols+i];
	C[(nrows-1)*ncols+i]=C[(nrows-2)*ncols+i];
	}
}