#include "NavierStokesCalc.h"
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
//Constructor
NavierStokesCalc::NavierStokesCalc(int inncols, int innrows, float indx, float indy, float inRe,  float inrho, float indynvis, bool inSMACbool, bool inSLICbool, bool inFLAPbool, bool inYVOFbool, float intmax, int inmaxNumIt, float inveleps, float inpeps, float ingx, float ingy, float inw, float intsf, float ingsf, bool inGSbool, bool inCGbool) : CFD(inncols,innrows,indx,indy,inRe,inrho,indynvis,inSMACbool,inSLICbool,inFLAPbool,inYVOFbool), maxNumIt(inmaxNumIt), tmax(intmax), veleps(inveleps), peps(inpeps), gx(ingx), gy(ingy), w(inw), tsf(intsf), gsf(ingsf), GSbool(inGSbool), CGbool(inCGbool)
{
}
NavierStokesCalc::~NavierStokesCalc()
{
}
//Methods
void NavierStokesCalc::RunNavierStokes()
{
DynamicAssignDTGAM();
ImposePeriodicVelocityBC();
//ImposevelocityBC();//ImposeSpecialVelocityBC();
//ImposepressureBC();
ImposePeriodicPressureBC();
FreeSurfaceConditions();
CalculatePreliminaryGuess();
//SOR:CalcRHS() or CG:
if(GSbool==true)
{
CalcRHS();
}
else if(CGbool==true)
{
InitialConditionsForCG2(); //or InitialConditionsForCG() ???
}
else
{
	//Print error code: No Pressure solution method chosen)
	cout<<"Error: No pressure solution method chosen"<<endl;
	exit(0);
}
bool pcheck=PressureResidualCheck();
int NoIt=0;
while((NoIt<maxNumIt) && (pcheck==true))
{
//SOR:SORpressure() or CG:CGpressure()
if(GSbool==true)
{
	SORpressure();
}
if(CGbool==true)
{
CGpressure();
}
//ImposepressureBC();
ImposePeriodicPressureBC();
pcheck=PressureResidualCheck();
NoIt+=1;
}
//ImposepressureBC();
ImposePeriodicPressureBC();
int mn=0;
}
void NavierStokesCalc::SetPreliminaryGuessValues()
{
	for(int j=0;j<nrows;j++)
	{
	for(int i=0;i<ncols;i++)
	{
		F.push_back(u[j*ncols+i]);
		G.push_back(v[j*ncols+i]);
		rhs.push_back(0.0);
		pd.push_back(0.0);
		r.push_back(0.0);
	}
	}
}
void NavierStokesCalc::CalculatePreliminaryGuess()
{
	for(int J=0;J<nrows;J++)
	{
	for(int I=0;I<ncols;I++)
	{
	 if((fl[J*ncols+I].C==true && fl[J*ncols+I].R==true && fl[J*ncols+I].L==true && fl[J*ncols+I].B==true && fl[J*ncols+I].T==false && fl[J*ncols+I].F==true)||(fl[J*ncols+I].C==true && fl[J*ncols+I].R==true && fl[J*ncols+I].L==false && fl[J*ncols+I].B==true && fl[J*ncols+I].T==false && fl[J*ncols+I].F==true))
		{
			//Boundary cell borders above or Boundary cell borders above and to the left
			Fcalc(I,J);
		}
		else if((fl[J*ncols+I].C==true && fl[J*ncols+I].R==false && fl[J*ncols+I].L==true && fl[J*ncols+I].B==true && fl[J*ncols+I].T==true && fl[J*ncols+I].F==true)||(fl[J*ncols+I].C==true && fl[J*ncols+I].R==false && fl[J*ncols+I].L==true && fl[J*ncols+I].B==false && fl[J*ncols+I].T==true && fl[J*ncols+I].F==true))
		{
			//Boundary Cell borders to the right or below and to the right
			Gcalc(I,J);
		}
		else if((fl[J*ncols+I].C==true && fl[J*ncols+I].R==true && fl[J*ncols+I].L==true && fl[J*ncols+I].B==true && fl[J*ncols+I].T==true && fl[J*ncols+I].F==true)||(fl[J*ncols+I].C==true && fl[J*ncols+I].R==true && fl[J*ncols+I].L==false && fl[J*ncols+I].B==true && fl[J*ncols+I].T==true && fl[J*ncols+I].F==true)||(fl[J*ncols+I].C==true && fl[J*ncols+I].R==true && fl[J*ncols+I].L==true && fl[J*ncols+I].B==false && fl[J*ncols+I].T==true && fl[J*ncols+I].F==true)||(fl[J*ncols+I].C==true && fl[J*ncols+I].R==true && fl[J*ncols+I].L==false && fl[J*ncols+I].B==false && fl[J*ncols+I].T==true && fl[J*ncols+I].F==true))
		{
			//Internal CV cell, Boundary to the left, boundary below or below and to the left
			Fcalc(I,J);
			Gcalc(I,J);
		}
	}
	}
}
void NavierStokesCalc::Fcalc(int &i, int &j)
{			
			float du2dxa=(1.0/dx)*((u[j*ncols+i]+u[j*ncols+i+1])/2.0)*((u[j*ncols+i]+u[j*ncols+i+1])/2.0)-(1.0/dx)*((u[j*ncols+i-1]+u[j*ncols+i])/2.0)*((u[j*ncols+i-1]+u[j*ncols+i])/2.0);
			float du2dxb=gam*(1.0/dx)*(((abs(u[j*ncols+i]+u[j*ncols+i+1]))/2.0)*((u[j*ncols+i]-u[j*ncols+i+1])/2.0)-(abs(u[j*ncols+i-1]+u[j*ncols+i])/2.0)*((u[j*ncols+i-1]-u[j*ncols+i])/2.0));
			float duvdya=(1.0/dy)*((u[j*ncols+i]+u[(j+1)*ncols+i])/2.0)*((v[j*ncols+i]+v[j*ncols+i+1])/2.0)-(1.0/dy)*((u[j*ncols+i]+u[(j-1)*ncols+i])/2.0)*((v[(j-1)*ncols+i]+v[(j-1)*ncols+i+1])/2.0);
			float duvdyb=gam*(1.0/dy)*((abs(v[j*ncols+i]+v[j*ncols+i+1])/2.0)*((u[j*ncols+i]-u[(j+1)*ncols+i])/2.0)-(abs(v[(j-1)*ncols+i]+v[(j-1)*ncols+i+1])/2.0)*((u[(j-1)*ncols+i]-u[j*ncols+i])/2.0));
			float CONX=du2dxa+du2dxb+duvdya+duvdyb;
			float VISX=(dynvis/rho)*((1.0/(dx*dx))*(u[j*ncols+i+1]-2.0*u[j*ncols+i]+u[j*ncols+i-1])+((1.0/(dy*dy))*(u[(j+1)*ncols+i]-2.0*u[j*ncols+i]+u[(j-1)*ncols+i])));
			F[j*ncols+i]=u[j*ncols+i]+dt*(VISX-CONX+gx);
}
void NavierStokesCalc::Gcalc(int &i, int &j)
{
			float dv2dya=(1.0/dy)*(((v[j*ncols+i]+v[(j+1)*ncols+i])/2.0)*((v[j*ncols+i]+v[(j+1)*ncols+i])/2.0)-((v[(j-1)*ncols+i]+v[j*ncols+i])/2.0)*((v[(j-1)*ncols+i]+v[j*ncols+i])/2.0));
			float dv2dyb=gam*(1.0/dy)*((abs(v[j*ncols+i]+v[(j+1)*ncols+i])/2.0)*((v[j*ncols+i]-v[(j+1)*ncols+i])/2.0)-(abs(v[(j-1)*ncols+i]+v[j*ncols+i])/2.0)*((v[(j-1)*ncols+i]-v[j*ncols+i])/2.0));
			float duvdxa=(1.0/dx)*(((u[j*ncols+i]+u[(j+1)*ncols+i])/2.0)*((v[j*ncols+i]+v[j*ncols+i+1])/2.0)-((u[j*ncols+i-1]+u[(j+1)*ncols+i-1])/2.0)*((v[j*ncols+i-1]+v[j*ncols+i])/2.0));
			float duvdxb=gam*(1.0/dy)*((abs(u[j*ncols+i]+u[(j+1)*ncols+i])/2.0)*((v[j*ncols+i]-v[j*ncols+i+1])/2.0)-(abs(u[j*ncols+i-1]+u[(j+1)*ncols+i-1])/2.0)*((v[j*ncols+i-1]-v[j*ncols+i])/2.0));
			float CONY=dv2dya+dv2dyb+duvdxa+duvdxb;
			float VISY=(dynvis/rho)*((1.0/(dx*dx))*(v[j*ncols+i+1]-2.0*v[j*ncols+i]+v[j*ncols+i-1])+((1.0/(dy*dy))*(v[(j+1)*ncols+i]-2.0*v[j*ncols+i]+v[(j-1)*ncols+i])));
			G[j*ncols+i]=v[j*ncols+i]+dt*(VISY-CONY+gy);
}
void NavierStokesCalc::ImposepressureBC()
{
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//CV borders above
			p[j*ncols+i]=p[(j+1)*ncols+i];
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//CV borders below
			p[j*ncols+i]=p[(j-1)*ncols+i];
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==false)
		{
			//CV borders left
			p[j*ncols+i]=p[j*ncols+i-1];
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==false)
		{
			//CV borders right
			p[j*ncols+i]=p[j*ncols+i+1];
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//CV borders above and to the right
			p[j*ncols+i]=(p[(j+1)*ncols+i]+p[j*ncols+i+1])/2.0;
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//CV borders above and to the left
			p[j*ncols+i]=(p[(j+1)*ncols+i]+p[j*ncols+i-1])/2.0;

		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//CV borders below and to the right
			p[j*ncols+i]=(p[(j-1)*ncols+i]+p[j*ncols+i+1])/2.0;

		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//CV borders below and to the left
			p[j*ncols+i]=(p[(j-1)*ncols+i]+p[j*ncols+i-1])/2.0;
		}
			//Boundary cell borders above
			//Boundary Cell borders to the right
			//Boundary cell borders above and to the left (equivalent to Boundary cell borders above)
			//Boundary cell borders below and to the right (equivalent to Boundary cell borders to the right)
			//Boundary cell borders above and to the right
			//:Do nothing
	}
}
}
void NavierStokesCalc::ImposeSpecialVelocityBC()
{
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//CV borders below, Impose Lid driving force of 1.0
			u[j*ncols+i]=2.0-u[(j-1)*ncols+i];
		}
	}
}
}
void NavierStokesCalc::ImposePeriodicVelocityBC()
{
	//Boundary Cells on far right equal boundary cells on far left, ie u[j*ncols+0]=u[j*ncols+ncols-2], with free slip at the vertical boundaries. 
	//Lower and upper Boundary Cells are still considered rigid, no slip.
	for(int j=0;j<nrows;j++)
	{
	u[j*ncols+0]=u[j*ncols+ncols-2];
	v[j*ncols+0]=v[j*ncols+ncols-2];
	v[j*ncols+ncols-1]=v[j*ncols+1];
	}
	for(int i=0;i<ncols;i++)
	{
	u[0*ncols+i]=-u[1*ncols+i];
	v[0*ncols+i]=0.0;
	u[(nrows-1)*ncols+i]=-u[(nrows-2)*ncols+i];
	v[(nrows-2)*ncols+i]=0.0;
	}		
}
void NavierStokesCalc::ImposePeriodicPressureBC()
{
	//Boundary Cells on far right equal boundary cells on far left, ie u[j*ncols+0]=u[j*ncols+ncols-2], with free slip at the vertical boundaries. 
	//Lower and upper Boundary Cells are still considered rigid, no slip.
	for(int j=0;j<nrows;j++)
	{
	p[j*ncols+0]=p[j*ncols+ncols-2];
	p[j*ncols+ncols-1]=p[j*ncols+1];
	}
	for(int i=0;i<ncols;i++)
	{
	p[0*ncols+i]=p[1*ncols+i];
	p[(nrows-1)*ncols+i]=p[(nrows-2)*ncols+i];
	}		
}
void NavierStokesCalc::ImposevelocityBC()
{
	//All Boundary Cells are considered rigid, no slip.
	for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//CV borders above
			u[j*ncols+i]=-u[(j+1)*ncols+i];
			v[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//CV borders below
			u[j*ncols+i]=-u[(j-1)*ncols+i];
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==false)
		{
			//CV borders left
			v[j*ncols+i]=-v[j*ncols+i-1];
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==false)
		{
			//CV borders right
			v[j*ncols+i]=-v[j*ncols+i+1];
			u[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//CV borders above and to the right
			u[j*ncols+i]=0.0;
			v[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//CV borders above and to the left
			v[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==false && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//CV borders below and to the right
			u[j*ncols+i]=0.0;
		}
			//CV borders below and to the left
			//Do nothing

		else if(fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==true && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//Boundary cell borders above
			v[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==true && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==true && fl[j*ncols+i].T==true)
		{
			//Boundary Cell borders to the right
			u[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false)
		{
			//Boundary cell borders above and to the left (equivalent to Boundary cell borders above)
			v[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==true && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//Boundary cell borders below and to the right (equivalent to Boundary cell borders to the right)
			u[j*ncols+i]=0.0;
		}
		else if(fl[j*ncols+i].C==true && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true)
		{
			//Boundary cell borders above and to the right
			u[j*ncols+i]=0.0;
			v[j*ncols+i]=0.0;
		}
	}
}
}
void NavierStokesCalc::CalcRHS()
{
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true)
		{
		rhs[j*ncols+i]=(rho/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy));
		}
	}
}
}
bool NavierStokesCalc::PressureResidualCheck()
{
//Pressure check calculates the residual (pressure gradients minus preliminary velocity gradients), calculates the L2norm and compares it to a tolerance.
//If the L2norm is greater than the tolerance the member function returns true, else it returns false. 
vector<float> err(ncols*nrows);
double sum=0;
double L2norm2=1.0/((ncols-2)*(nrows-2));
for(int j=0;j<nrows;j++)
{
for(int i=0;i<ncols;i++)
{
if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true)
		{
		err[j*ncols+i]=((p[j*ncols+i+1]-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))+((p[(j+1)*ncols+i]-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy))-rhs[j*ncols+i];
		//r3[j*ncols+i]=abs(r[j*ncols+i]);
		sum+=err[j*ncols+i]*err[j*ncols+i];
		}
}
}
double L2norm=sqrt(L2norm2*sum);
if(L2norm>peps)
{
return true;
}
else
{
return false;
}
}

void NavierStokesCalc::SORpressure()
{

for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true)
		{
		p[j*ncols+i]=(1.0-w)*p[j*ncols+i]+(w/((2.0/(dx*dx))+(2.0/(dy*dy))))*(((p[j*ncols+i+1]+p[j*ncols+i-1])/(dx*dx))+((p[(j+1)*ncols+i]+p[(j-1)*ncols+i])/(dy*dy))-rhs[j*ncols+i]);
		}
	}
}
}
void NavierStokesCalc::InitialConditionsForCG()
{
//Calculate adapted rhs (where non changeable pressures, pressures not being solved for in the system of linear equations (in ghost cells and empty cells), are subtracted from the rhs
//Calculate initial values for pd and r
rsum=0.0;
sumpd=0.0;
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//Fluid on all sides
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy));
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i+1]-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))+((p[(j+1)*ncols+i]-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));

		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid on the right	
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dx*dx))*p[j*ncols+i+1];
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i-1]-2.0*p[j*ncols+i])/(dx*dx))-((p[(j+1)*ncols+i]-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));

		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid on the left
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dx*dx))*p[j*ncols+i-1];
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i+1]-2.0*p[j*ncols+i])/(dx*dx))-((p[(j+1)*ncols+i]-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dy*dy))*p[(j+1)*ncols+i];
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i+1]-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))+((-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));

		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below	
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dy*dy))*p[(j-1)*ncols+i];
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i+1]-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))+((p[(j+1)*ncols+i]-2.0*p[j*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above and to the right	
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dy*dy))*p[(j+1)*ncols+i]-(1.0/(dx*dx))*p[j*ncols+i+1];
		r[j*ncols+i]=rhs[j*ncols+i]-((-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))-((-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below and on the right	
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dy*dy))*p[(j-1)*ncols+i]-(1.0/(dx*dx))*p[j*ncols+i+1];
		r[j*ncols+i]=rhs[j*ncols+i]-((-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))-((p[(j+1)*ncols+i]-2.0*p[j*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above and on the left	
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dy*dy))*p[(j+1)*ncols+i]-(1.0/(dx*dx))*p[j*ncols+i-1];
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i+1]-2.0*p[j*ncols+i])/(dx*dx))-((-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+1]*r[j*ncols+i]);
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below and on the left	
		rhs[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-(1.0/(dy*dy))*p[(j-1)*ncols+i]-(1.0/(dx*dx))*p[j*ncols+i-1];
		r[j*ncols+i]=rhs[j*ncols+i]-((p[j*ncols+i+1]-2.0*p[j*ncols+i])/(dx*dx))-((p[(j+1)*ncols+i]-2.0*p[j*ncols+i])/(dy*dy));
		pd[j*ncols+i]=r[j*ncols+i];
		rsum+=(r[j*ncols+i]*r[j*ncols+i]);
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
	}
}
}
void NavierStokesCalc::InitialConditionsForCG2()
{
//Calculate adapted rhs (where non changeable pressures, pressures not being solved for in the system of linear equations (in ghost cells and empty cells), are subtracted from the rhs
//Calculate initial values for pd and r
rsum=0.0;
sumpd=0.0;
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true)
		{
			r[j*ncols+i]=(1.0/dt)*(((F[j*ncols+i]-F[j*ncols+i-1])/dx)+((G[j*ncols+i]-G[(j-1)*ncols+i])/dy))-((p[j*ncols+i+1]-2.0*p[j*ncols+i]+p[j*ncols+i-1])/(dx*dx))+((p[(j+1)*ncols+i]-2.0*p[j*ncols+i]+p[(j-1)*ncols+i])/(dy*dy));
			pd[j*ncols+i]=r[j*ncols+i];
			rsum+=(r[j*ncols+i]*r[j*ncols+i]);
			if(fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
			{
			//Fluid on all sides
			sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
			{
			//No Fluid on the right	
			sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
			{
			//No Fluid on the left
			sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
			{
			//No Fluid above
			sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
			{
			//No Fluid below	
			sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
			{
			//No Fluid above and to the right	
			sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
			{
			//No Fluid below and on the right	
			sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
			{
			//No Fluid above and on the left	
			sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
			else if(fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
			{
			//No Fluid below and on the left	
			sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
			}
		}
	}
}
}
void NavierStokesCalc::CGpressure()
{

alpha=rsum/sumpd;
newrsum=0;
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//Fluid on all sides
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((pd[j*ncols+i+1]-2.0*pd[j*ncols+i]+pd[j*ncols+i-1])/(dx*dx))+((pd[(j+1)*ncols+i]-2.0*pd[j*ncols+i]+pd[(j-1)*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid on the right	
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((-2.0*pd[j*ncols+i]+pd[j*ncols+i-1])/(dx*dx))+((pd[(j+1)*ncols+i]-2.0*pd[j*ncols+i]+pd[(j-1)*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid on the left
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((pd[j*ncols+i+1]-2.0*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]-2.0*pd[j*ncols+i]+pd[(j-1)*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((pd[j*ncols+i+1]-2.0*pd[j*ncols+i]+pd[j*ncols+i-1])/(dx*dx))+((-2.0*pd[j*ncols+i]+pd[(j-1)*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below	
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((pd[j*ncols+i+1]-2.0*pd[j*ncols+i]+pd[j*ncols+i-1])/(dx*dx))+((pd[(j+1)*ncols+i]-2.0*pd[j*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above and to the right	
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((-2.0*pd[j*ncols+i]+pd[j*ncols+i-1])/(dx*dx))+((-2.0*pd[j*ncols+i]+pd[(j-1)*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below and on the right	
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((-2.0*pd[j*ncols+i]+pd[j*ncols+i-1])/(dx*dx))+((pd[(j+1)*ncols+i]-2.0*pd[j*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above and on the left	
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((pd[j*ncols+i+1]-2.0*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]+pd[(j-1)*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below and on the left	
		p[j*ncols+i]+=alpha*pd[j*ncols+i];
		r[j*ncols+i]-=alpha*((pd[j*ncols+i+1]-2.0*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]-2.0*pd[j*ncols+i])/(dy*dy));
		newrsum+=(r[j*ncols+i]*r[j*ncols+i]);
		}
	}
}
beta=newrsum/rsum;
rsum=newrsum;
sumpd=0;
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true)
		{
		//Update direction vector pd
		pd[j*ncols+i]=r[j*ncols+i]+beta*pd[j*ncols+i];
		}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//Fluid on all sides
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid on the right	
		sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid on the left
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below	
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above and to the right	
		sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==false && fl[j*ncols+i-1].F==true && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below and on the right	
		sumpd+=((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[j*ncols+i-1]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==false && fl[(j-1)*ncols+i].F==true)
		{
		//No Fluid above and on the left	
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((-2.0*pd[j*ncols+i]*pd[j*ncols+i]+pd[(j-1)*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
		else if(fl[j*ncols+i].C==true  && fl[j*ncols+i].F==true && fl[j*ncols+i+1].F==true && fl[j*ncols+i-1].F==false && fl[(j+1)*ncols+i].F==true && fl[(j-1)*ncols+i].F==false)
		{
		//No Fluid below and on the left	
		sumpd+=((pd[j*ncols+i+1]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dx*dx))+((pd[(j+1)*ncols+i]*pd[j*ncols+i]-2.0*pd[j*ncols+i]*pd[j*ncols+i])/(dy*dy));
		}
	}
}
}
bool NavierStokesCalc::UpdateVelocity()
{
double L2norm2=1.0/((ncols-2)*(nrows-2));
double sumx=0;
double sumy=0;
for(int j=0;j<nrows;j++)
{
		for(int i=0;i<ncols;i++)
		{
			if((fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==true && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false && fl[j*ncols+i].F==true)||(fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==false && fl[j*ncols+i].F==true))
			{
			F[j*ncols+i]-=(dt/dx)*(1.0/rho)*(p[j*ncols+i+1]-p[j*ncols+i]);
			sumx+=(F[j*ncols+i]-u[j*ncols+i])*(F[j*ncols+i]-u[j*ncols+i]);
			u[j*ncols+i]=F[j*ncols+i];
			}
			else if((fl[j*ncols+i].C==true && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==true && fl[j*ncols+i].T==true && fl[j*ncols+i].F==true)||(fl[j*ncols+i].C==true && fl[j*ncols+i].R==false && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true && fl[j*ncols+i].F==true))
			{
			G[j*ncols+i]-=(dt/dy)*(1.0/rho)*(p[(j+1)*ncols+i]-p[j*ncols+i]);
			sumy+=(G[j*ncols+i]-v[j*ncols+i])*(G[j*ncols+i]-v[j*ncols+i]);
			v[j*ncols+i]=G[j*ncols+i];
			}
			else if((fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==true && fl[j*ncols+i].B==true && fl[j*ncols+i].T==true && fl[j*ncols+i].F==true)||(fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==true && fl[j*ncols+i].T==true && fl[j*ncols+i].F==true)||(fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==true && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true && fl[j*ncols+i].F==true)||(fl[j*ncols+i].C==true && fl[j*ncols+i].R==true && fl[j*ncols+i].L==false && fl[j*ncols+i].B==false && fl[j*ncols+i].T==true && fl[j*ncols+i].F==true))
			{
			F[j*ncols+i]-=(dt/dx)*(1.0/rho)*(p[j*ncols+i+1]-p[j*ncols+i]);
			sumx+=(F[j*ncols+i]-u[j*ncols+i])*(F[j*ncols+i]-u[j*ncols+i]);
			u[j*ncols+i]=F[j*ncols+i];
			G[j*ncols+i]-=(dt/dy)*(1.0/rho)*(p[(j+1)*ncols+i]-p[j*ncols+i]);
			sumy+=(G[j*ncols+i]-v[j*ncols+i])*(G[j*ncols+i]-v[j*ncols+i]);
			v[j*ncols+i]=G[j*ncols+i];
			}
		}
}
double L2normx=sqrt(L2norm2*sumx);
double L2normy=sqrt(L2norm2*sumy);
if(L2normx>veleps || L2normy>veleps)
{
return true;
}
else
{
return false;
}
}
void NavierStokesCalc::DynamicAssignDTGAM()
{
vector<float> stabt(4);
vector<float> stabg(3);
vector<float>::iterator it;
vector<float> ua(ncols*nrows);
vector<float> va(ncols*nrows);
for(int j=1;j<=(nrows-2);j++)
{
	for(int i=1;i<=(ncols-2);i++)
	{
		ua[j*ncols+i]=abs(u[j*ncols+i]);
		va[j*ncols+i]=abs(v[j*ncols+i]);
	}
}
it=max_element(ua.begin(),ua.end());
stabt[0]=dx/(*it);
it=max_element(va.begin(),va.end());
stabt[1]=dy/(*it);
stabt[2]=(Re/2.0)*(1.0/((1.0/(dx*dx))+(1.0/(dy*dy))));
stabt[3]=0.1;
it=min_element(stabt.begin(),stabt.begin()+4);
dt=(*it)*tsf;
it=max_element(ua.begin(),ua.end());
stabg[0]=(*it)*dt/dx;
it=max_element(va.begin(),va.end());
stabg[1]=(*it)*dt/dy;
it=max_element(stabg.begin(),stabg.begin()+3);
stabg[2]=0.1;
gam=(*it)*gsf;
}
