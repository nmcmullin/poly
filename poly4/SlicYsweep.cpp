#include "FreeSurfaceCalc.h"
#include "CFD.h"
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

void FreeSurfaceCalc::SLICYsweep(int &I, int &J)
{
float width;
C[J*ncols+I];
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(C[J*ncols+I]<(1.0-voftol) && (1.0-C[J*ncols+I])<(1.0-voftol))//i.e if the cell contains an interface
	{
		if((C[(J+1)*ncols+I]>voftol && (1.0-C[(J+1)*ncols+I])>voftol && C[(J-1)*ncols+I]>voftol && (1-C[(J-1)*ncols+I])>voftol ) || (C[(J+1)*ncols+I]<=voftol && (1-C[(J+1)*ncols+I])<=voftol && C[(J-1)*ncols+I]<=voftol && (1-C[(J-1)*ncols+I])<=voftol) || (C[(J+1)*ncols+I]<=voftol && (1-C[(J+1)*ncols+I])<=voftol && C[(J-1)*ncols+I]>voftol && (1-C[(J-1)*ncols+I])>voftol) || (C[(J+1)*ncols+I]>voftol && (1-C[(J+1)*ncols+I])>voftol && C[(J-1)*ncols+I]<=voftol && (1-C[(J-1)*ncols+I])<=voftol))
		{
		//Type I, vertical interface  equivalent to the type 1 horizontal interface in x-sweep
		width=C[J*ncols+I];
			if(v[J*ncols+I]>0.0)
			{
				Ft[J*ncols+I]=v[J*ncols+I]*dt*width*dx;
			}
			if(v[(J-1)*ncols+I]<0.0)
			{			
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*width*dx;
			}
		}
		else if((C[(J-1)*ncols+I]>voftol && C[(J+1)*ncols+I]<=voftol)||(C[(J-1)*ncols+I]<=voftol && C[(J+1)*ncols+I]>voftol)||((1.0-C[(J-1)*ncols+I])>voftol && (1.0-C[(J+1)*ncols+I])<=voftol)||((1.0-C[(J-1)*ncols+I])<=voftol && (1.0-C[(J+1)*ncols+I])>voftol))
		{
		//Type II, horizontal interface 
		width=C[J*ncols+I];
		if(C[(J-1)*ncols+i]>C[(J+1)*ncols+i])
		{
		//Highlighted Fluid below
			//equivalent st=0.0	sb=1.0;
			if(v[J*ncols+I]>0.0)
			{
				if(v[J*ncols+I]*dt>(1-width)*dy)
				{
				Ft[J*ncols+I]=(v[J*ncols+I]*dt-(1.0-width)*dy)*dx;
				}
				else
				{
				Ft[J*ncols+I]=0.0;
				}
			}
			if(v[(J-1)*ncols+I]<0.0)
			{
				if(abs(v[(J-1)*ncols+I])*dt<width*dy)
				{
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx;
				}
				else
				{
				Fb[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
		}
		else
		{
		// Highlighted Fluid Above
			//equivalent st=1.0	sb=0.0;
			if(v[J*ncols+I]>0.0)
			{
				if(v[J*ncols+I]*dt<width*dy)
				{
				Ft[J*ncols+I]=v[J*ncols+I]*dt*dx;
				}
				else
				{
				Ft[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
			if(v[(J-1)*ncols+I]<0.0)
			{
				if(abs(v[(J-1)*ncols+I])*dt>(1.0-width)*dy)
				{
				Fb[J*ncols+I]=(abs(v[(J-1)*ncols+I])*dt-(1.0-width)*dy)*dx;
				}
				else
				{
				Fb[J*ncols+I]=0.0;
				}
			}
		}
		}
		else if((C[(J+1)*ncols+I]<=voftol && C[(J-1)*ncols+I]<=voftol && (1.0-C[(J+1)*ncols+I])>voftol && (1.0-C[(J-1)*ncols+I])>voftol)||(C[(J+1)*ncols+I]>voftol && C[(J-1)*ncols+I]>voftol && (1.0-C[(J+1)*ncols+I])<=voftol && (1.0-C[(J-1)*ncols+I])<=voftol))
		{
		//Type III, two horizontal interfaces (sandwiching a fluid between them). The sandwiched fluid is centered 
		width=C[J*ncols+I];
		if(C[(J+1)*ncols+I]>0.0 && C[(J-1)*ncols+I]>0.0)
		{
		//equivalent st=1.0	sb=1.0;
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<0.5*width*dy)
			{
				Ft[J*ncols+I]=v[J*ncols+i]*dt*dx;
			}
			else if(v[J*ncols+I]*dt<(0.5*width*dy+(1.0-width)*dy))
			{
				Ft[J*ncols+I]=0.5*width*dx*dy;
			}
			else
			{
				Ft[J*ncols+I]=0.5*width*dx*dy+(v[J*ncols+I]*dt-(0.5*width*dy+(1.0-width)*dy))*dx;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<0.5*width*dy)
			{
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx;
			}
			else if(abs(v[(J-1)*ncols+I])*dt<(0.5*width*dy+(1.0-width)*dy))
			{
				Fb[J*ncols+I]=0.5*width*dx*dy;
			}
			else
			{
				Fb[J*ncols+I]=0.5*width*dx*dy+(abs(v[(J-1)*ncols+I])*dt-(0.5*width*dy+(1.0-width)*dy))*dx;
			}
		}	
		}
		else
		{
		//equivalent of st=0.0;sb=0.0;
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<0.5*(1.0-width)*dy)
			{
				Ft[J*ncols+I]=0.0;
			}
			else if(v[J*ncols+I]*dt<(0.5*(1.0-width)*dy+width*dy))
			{
				Ft[J*ncols+I]=(v[J*ncols+I]*dt-0.5*(1.0-width)*dy)*dx;
			}
			else
			{
				Ft[J*ncols+I]=width*dx*dy;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<0.5*(1.0-width)*dy)
			{
				Fb[J*ncols+I]=0.0;
			}
			else if(abs(v[(J-1)*ncols+I])*dt<(0.5*(1.0-width)*dy+width*dy))
			{
				Fb[J*ncols+I]=(abs(v[(J-1)*ncols+I])*dt-0.5*(1.0-width)*dy)*dx;
			}
			else
			{
				Fb[J*ncols+I]=width*dx*dy;
			}
		}
		}
		}
	}
	}
}
}