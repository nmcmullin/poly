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

void FreeSurfaceCalc::SLICXsweep(int &I, int &J)
{
float width;
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(C[J*ncols+I]<(1.0-voftol) && C[J*ncols+I]>voftol)//(1.0-C[J*ncols+I])<(1.0-voftol))//i.e if the cell contains an interface
	{
		if((C[J*ncols+I+1]>voftol && (1.0-C[J*ncols+I+1])>voftol && C[J*ncols+I-1]>voftol && (1-C[J*ncols+I-1])>voftol ) || (C[J*ncols+I+1]<=voftol && (1-C[J*ncols+I+1])<=voftol && C[J*ncols+I-1]<=voftol && (1-C[J*ncols+I-1])<=voftol) || (C[J*ncols+I+1]<=voftol && (1-C[J*ncols+I+1])<=voftol && C[J*ncols+I-1]>voftol && (1-C[J*ncols+I-1])>voftol) || (C[J*ncols+I+1]>voftol && (1-C[J*ncols+I+1])>voftol && C[J*ncols+I-1]<=voftol && (1-C[J*ncols+I-1])<=voftol))
		{
		//Type I, horizontal interface with the interface at a height of the VF*dy of the lower fluid (the one with VF=1 below) 
			width=C[J*ncols+I];
			if(u[J*ncols+I]>0.0)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*width*dy;
			}
			if(u[J*ncols+I-1]<0.0)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*width*dy;
			}	
		}
		else if((C[J*ncols+I-1]>voftol && C[J*ncols+I+1]<=voftol)||(C[J*ncols+I-1]<=voftol && C[J*ncols+I+1]>voftol)||((1.0-C[J*ncols+I-1])>voftol && (1.0-C[J*ncols+I+1])<=voftol)||((1.0-C[J*ncols+I-1])<=voftol && (1.0-C[J*ncols+I+1])>voftol))
		{
		//Type II, vertical interface located at VF*dx away from the cells left boundary where VF is the left hand fluid (the fluid with VF==1 in the cell to the left
		width=C[J*ncols+I];
		if(C[J*ncols+i-1]>C[J*ncols+i+1])
		{
			//equivalent sr=0.0 sl=1.0;
			if(u[J*ncols+I]>0)
			{
				if(u[J*ncols+I]*dt>(1.0-width)*dx)
				{
				Fr[J*ncols+I]=(u[J*ncols+I]*dt-(1.0-width)*dx)*dy;
				}
				else
				{
				Fr[J*ncols+I]=0.0;
				}
			}
			if(u[J*ncols+I-1]<0.0)
			{			
				if(abs(u[J*ncols+I-1])*dt<width*dx)//i-->i-1********
				{
				Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
				}
				else
				{
				Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
		}
		else
		{
		//equivalent sr=1.0 sl=0.0;
			if(u[J*ncols+I]>0)
			{
				if(u[J*ncols+I]*dt<width*dx)
				{
				Fr[J*ncols+I]=u[J*ncols+I]*dt*dy;
				}
				else
				{
				Fr[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
			if(u[J*ncols+I-1]<0.0)
			{			
				if(abs(u[J*ncols+I-1])*dt>(1.0-width)*dx)//i-->i-1*****
				{
				Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-(1-width)*dx)*dy;
				}
				else
				{
				Fl[J*ncols+I]=0.0;
				}
			}
		}
		}
		else if((C[J*ncols+I+1]<=voftol && C[J*ncols+I-1]<=voftol && (1.0-C[J*ncols+I+1])>voftol && (1.0-C[J*ncols+I+1])>voftol)||(C[J*ncols+I+1]>voftol && C[J*ncols+I-1]>voftol && (1.0-C[J*ncols+I+1])<=voftol && (1.0-C[J*ncols+I-1])<=voftol))
		{
		//Type III, two vertical interfaces (sandwiching a fluid between them). Presumably the sandwiched fluid is centered however I am uncertain about this
		width=C[J*ncols+I];
		if(C[J*ncols+I+1]>0 && C[J*ncols+I-1]>0)
		{
		//equivalent sr=1.0 sl=1.0;
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<0.5*width*dx)
			{
				Fr[J*ncols+I]=u[J*ncols+i]*dt*dy;
			}
			else if(u[J*ncols+I]*dt<(0.5*width*dx+(1.0-width)*dx))
			{
				Fr[J*ncols+I]=0.5*width*dx*dy;
			}
			else
			{
				Fr[J*ncols+I]=0.5*width*dx*dy+(u[J*ncols+I]*dt-(0.5*width*dx+(1.0-width)*dx))*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<0.5*width*dx)
			{
				Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
			}
			else if(abs(u[J*ncols+I-1])*dt<(0.5*width*dx+(1.0-width)*dx))
			{
				Fl[J*ncols+I]=0.5*width*dx*dy;
			}
			else
			{
				Fl[J*ncols+I]=0.5*width*dx*dy+(abs(u[J*ncols+I-1])*dt-(0.5*width*dx+(1.0-width)*dx))*dy;
			}
		}	
		}
		else
		{
		//equivalent sr=0.0 sl=0.0;
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<0.5*(1.0-width)*dx)
			{
				Fr[J*ncols+I]=0.0;
			}
			else if(u[J*ncols+I]*dt<(0.5*(1.0-width)*dx+width*dx))
			{
				Fr[J*ncols+I]=(u[J*ncols+I]*dt-0.5*(1.0-width)*dx)*dy;
			}
			else
			{
				Fr[J*ncols+I]=width*dx*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<0.5*(1.0-width)*dx)
			{
				Fl[J*ncols+I]=0.0;
			}
			else if(abs(u[J*ncols+I-1])*dt<(0.5*(1.0-width)*dx+width*dx))
			{
				Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-0.5*(1.0-width)*dx)*dy;
			}
			else
			{
				Fl[J*ncols+I]=width*dx*dy;
			}
		}
		}
		}
	}
	}
}
}