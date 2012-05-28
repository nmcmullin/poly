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

void FreeSurfaceCalc::FluxCalcXsweepC(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL)
{
if(ALPHA<(PI/4))
{
	if(C[J*ncols+I]<=(0.5*tan(ALPHA)))
	{
	//case I
	if(u[J*ncols+I]>0.0)
	{
		if(u[J*ncols+I]*dt<=(1.0-ST)*dx)
		{
		Fr[J*ncols+I]=0.0;
		}
		else
		{
		Fr[J*ncols+I]=0.5*(u[J*ncols+I]*dt-(1.0-ST)*dx)*(u[J*ncols+I]*dt-(1.0-ST)*dx)*tan(BETA);
		}
	}
	if(u[J*ncols+I-1]<0.0)
	{
		if(abs(u[J*ncols+I-1])*dt>=ST*dx)
		{
		Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
		}
		else
		{
		Fl[J*ncols+I]=0.5*abs(u[J*ncols+I-1])*dt*(2.0-abs(u[J*ncols+I-1])*dt/(ST*dx))*SL*dy;
		}
	}
	}
	else if(C[J*ncols+I]<=(1.0-0.5*tan(ALPHA)))
	{
	//Case II
	if(u[J*ncols+I]>0.0)
	{
	Fr[J*ncols+I]=u[J*ncols+I]*dt*(SR*dy-0.5*u[J*ncols+I]*dt*tan(BETA));
	}
	if(u[J*ncols+I-1]<0.0)
	{
	Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*(SL*dy+0.5*abs(u[J*ncols+I-1])*dt*tan(BETA));
	}
	}
	else
	{
	//Case IV
	if(u[J*ncols+I]>0.0)
	{
		if(u[J*ncols+I]*dt<=(1.0-SB)*dx)
		{
		Fr[J*ncols+I]=u[J*ncols+I]*dt*dy-(1.0-C[J*ncols+I])*dx*dy;
		}
		else
		{
		Fr[J*ncols+I]=u[J*ncols+I]*dt*(SR*dy+0.5*u[J*ncols+I]*dt*tan(BETA));
		}
	}
	if(u[J*ncols+I-1]<0.0)
	{
		if(abs(u[J*ncols+I-1])*dt<=SB*dx)
		{
		Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
		}
		else
		{
		Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy-0.5*tan(BETA)*(abs(u[J*ncols+I-1])*dt-SB*dx)*(abs(u[J*ncols+I-1])*dt-SB*dx);
		}
	}
	}
}
else
{
	if(C[J*ncols+I]<=(0.5*(1.0/tan(ALPHA))))
	{
	//case I
	if(u[J*ncols+I]>0.0)
	{
		if(u[J*ncols+I]*dt<=(1.0-ST)*dx)
		{
		Fr[J*ncols+I]=0.0;
		}
		else
		{
		Fr[J*ncols+I]=0.5*(u[J*ncols+I]*dt-(1.0-ST)*dx)*(u[J*ncols+I]*dt-(1.0-ST)*dx)*tan(BETA);
		}
	}
	if(u[J*ncols+I-1]<0.0)
	{
		if(abs(u[J*ncols+I-1])*dt>=ST*dx)
		{
		Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
		}
		else
		{
		Fl[J*ncols+I]=0.5*abs(u[J*ncols+I-1])*dt*(2.0-abs(u[J*ncols+I-1])*dt/(ST*dx))*SL*dy;
		}
	}
	}
	else if(C[J*ncols+I]<=(1.0-0.5*(1.0/tan(ALPHA))))
	{
	//Case III
	if(u[J*ncols+I]>0.0)
	{
		if(u[J*ncols+I]*dt<=ST*dx)
		{
		Fr[J*ncols+I]=0.0;
		}
		else if(u[J*ncols+I]*dt<=SB*dx)
		{
		Fr[J*ncols+I]=0.5*(u[J*ncols+I]*dt-(1.0-ST)*dx)*(u[J*ncols+I]*dt-(1.0-ST)*dx)*tan(BETA);
		}
		else
		{
		Fr[J*ncols+I]=u[J*ncols+I]*dt*dy-(1.0-C[J*ncols+I])*dx*dy;
		}
	}
	if(u[J*ncols+I-1]<0.0)
	{
		if(abs(u[J*ncols+I-1])*dt<=SB*dx)
		{
		Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
		}
		else if(u[J*ncols+I]*dt<=ST*dx)
		{
		Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy-0.5*(abs(u[J*ncols+I-1])*dt-SB*dx)*(abs(u[J*ncols+I-1])*dt-SB*dx)*tan(BETA);//*********
		}
		else
		{
		Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
		}
	}
	}
	else
	{
	//Case IV
	if(u[J*ncols+I]>0.0)
	{
		if(u[J*ncols+I]*dt<=(1.0-SB)*dx)
		{
		Fr[J*ncols+I]=u[J*ncols+I]*dt*dy-(1.0-C[J*ncols+I])*dx*dy;
		}
		else
		{
		Fr[J*ncols+I]=u[J*ncols+I]*dt*(SR*dy+0.5*u[J*ncols+I]*dt*tan(BETA));
		}
	}
	if(u[J*ncols+I-1]<0.0)
	{
		if(abs(u[J*ncols+I-1])*dt<=SB*dx)
		{
		Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
		}
		else
		{
		Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy-0.5*tan(BETA)*(abs(u[J*ncols+I-1])*dt-SB*dx)*(abs(u[J*ncols+I-1])*dt-SB*dx);
		}
	}
	}
}
}