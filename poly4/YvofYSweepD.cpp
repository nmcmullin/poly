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
void FreeSurfaceCalc::FluxCalcYsweepDIIrevIII(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL)
{
if(ALPHA<(PI/4))
{
	if(C[J*ncols+I]<=(0.5*tan(ALPHA)))
	{
	//case I
	if(v[J*ncols+I]>0.0)
	{
		if(v[J*ncols+I]*dt>=SR*dy)
		{
		Ft[J*ncols+I]=C[J*ncols+I]*dx*dy;
		}
		else
		{
		Ft[J*ncols+I]=0.5*(v[J*ncols+I]*dt*(2.0-v[J*ncols+I]*dt/(SR*dy)))*ST*dx;
		}
	}
	if(v[(J-1)*ncols+I]<0.0)
	{
		if(abs(v[(J-1)*ncols+I])*dt<=(1.0-SR)*dy)
		{
		Fb[J*ncols+I]=0.0;
		}
		else
		{
		Fb[J*ncols+I]=0.5*(abs(v[(J-1)*ncols+I])*dt-(1.0-SR)*dy)*(abs(v[(J-1)*ncols+I])*dt-(1.0-SR)*dy)*(1.0/tan(BETA));
		}		
	}
	}
	else if(C[J*ncols+I]<=(1.0-0.5*tan(ALPHA)))
	{
	//Case III
	if(v[J*ncols+I]>0.0)
	{
	Ft[J*ncols+I]=v[J*ncols+I]*dt*(ST*dx-0.5*v[J*ncols+I]*dt*(1.0/tan(BETA)));
	}
	if(v[(J-1)*ncols+I]<0.0)
	{
	Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*(SB*dx+0.5*abs(v[(J-1)*ncols+I])*dt*(1.0/tan(BETA)));
	}
	}
	else
	{
	//Case IV
	if(v[J*ncols+I]>0.0)
	{
		if(v[J*ncols+I]*dt<=SL*dy)
		{
		Ft[J*ncols+I]=v[J*ncols+I]*dt*dx;
		}
		else
		{
		Ft[J*ncols+I]=v[J*ncols+I]*dt*dx-(0.5*(v[J*ncols+I]*dt-SL*dy)*(v[J*ncols+I]*dt-SL*dy)*(1.0/tan(BETA)));
		}
	}
	if(v[(J-1)*ncols+I]<0.0)
	{
		if(abs(v[(J-1)*ncols+I])*dt>=(1.0-SL)*dy)
		{
		Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx-(1.0-C[J*ncols+I])*dx*dy;
		}
		else
		{
		Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*(SB*dx+0.5*abs(v[(J-1)*ncols+I])*dt*(1.0/tan(BETA)));
		}
	}
	}
}
else
{
	if(C[J*ncols+I]<=(0.5*(1.0/tan(ALPHA))))
	{
	//case I
	if(v[J*ncols+I]>0.0)
	{
		if(v[J*ncols+I]*dt>=SR*dy)
		{
		Ft[J*ncols+I]=C[J*ncols+I]*dx*dy;
		}
		else
		{
		Ft[J*ncols+I]=0.5*(v[J*ncols+I]*dt*(2.0-v[J*ncols+I]*dt/(SR*dy)))*ST*dx;
		}
	}
	if(v[(J-1)*ncols+I]<0.0)
	{
		if(abs(v[(J-1)*ncols+I])*dt<=(1.0-SR)*dy)
		{
		Fb[J*ncols+I]=0.0;
		}
		else
		{
		Fb[J*ncols+I]=0.5*(abs(v[(J-1)*ncols+I])*dt-(1.0-SR)*dy)*(abs(v[(J-1)*ncols+I])*dt-(1.0-SR)*dy)*(1.0/tan(BETA));
		}		
	}	
	}
	else if(C[J*ncols+I]<=(1.0-0.5*(1.0/tan(ALPHA))))
	{
		//Case II
	if(v[J*ncols+I]>0.0)
	{
		if(v[J*ncols+I]*dt<=SL*dy)
		{
		Ft[J*ncols+I]=v[J*ncols+I]*dt*dx;
		}
		else if(v[J*ncols+I]*dt<=SR*dy)
		{
		Ft[J*ncols+I]=v[J*ncols+I]*dt*dx-0.5*(v[J*ncols+I]*dt-SL*dy)*(v[J*ncols+I]*dt-SL*dy)*(1.0/tan(BETA));
		}
		else
		{
		Ft[J*ncols+I]=C[J*ncols+I]*dx*dy;
		}
	}
	if(v[(J-1)*ncols+I]<0.0)
	{
		if(abs(v[(J-1)*ncols+I])*dt<=(1.0-SR)*dy)
		{
		Fb[J*ncols+I]=0.0;
		}
		else if(abs(v[(J-1)*ncols+I])*dt<=(1.0-SL)*dy)
		{
		Fb[J*ncols+I]=0.5*(abs(v[(J-1)*ncols+I])*dt-(1.0-SR)*dy)*(abs(v[(J-1)*ncols+I])*dt-(1.0-SR)*dy)*(1.0/tan(BETA));
		}
		else
		{
		Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx-(1.0-C[J*ncols+I])*dx*dy;
		}
	}
	}
	else
	{
	//Case IV
	if(v[J*ncols+I]>0.0)
	{
		if(v[J*ncols+I]*dt<=SL*dy)
		{
		Ft[J*ncols+I]=v[J*ncols+I]*dt*dx;
		}
		else
		{
		Ft[J*ncols+I]=v[J*ncols+I]*dt*dx-(0.5*(v[J*ncols+I]*dt-SL*dy)*(v[J*ncols+I]*dt-SL*dy)*(1.0/tan(BETA)));
		}
	}
	if(v[(J-1)*ncols+I]<0.0)
	{
		if(abs(v[(J-1)*ncols+I])*dt>=(1.0-SL)*dy)
		{
		Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx-(1.0-C[J*ncols+I])*dx*dy;
		}
		else
		{
		Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*(SB*dx+0.5*abs(v[(J-1)*ncols+I])*dt*(1.0/tan(BETA)));
		}
	}
	}
}
}