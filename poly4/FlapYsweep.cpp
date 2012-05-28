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

void FreeSurfaceCalc::FLAPYsweep(int &I, int &J)
{
	float sr,sl,st,sb;
//Assuming the cell geometry and location determinations remain the same while calculating the y direction fluxes. i.e PLIC
	if((C[(J-1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]>voftol) && (C[(J+1)*ncols+I]<(1.0-voftol) && C[(J+1)*ncols+I]>voftol))
	{
		//Vertical Interface
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		if(C[(J-1)*ncols+I]>C[(J+1)*ncols+I])
		{	
			//Vert B on L
			sr=0.0;
			sl=1.0;
			if(v[J*ncols+I]>0.0)
			{
				Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
			}
			if(v[(J-1)*ncols+I]<0.0)
			{
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*sb*dx;
			}
		}
		else
		{
			//Vert B on R
			sr=1.0;
			sl=0.0;
			if(v[J*ncols+I]>0.0)
			{
				Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
			}
			if(v[(J-1)*ncols+I]<0.0)
			{
				Fb[J*ncols+I]=abs(v[J*ncols+I])*dt*sb*dx;
			}
		}
	}
	else if(C[(J-1)*ncols+I]>voftol  && C[(J+1)*ncols+I]<=voftol)
	{
	if((C[J*ncols+I+1]>voftol && C[J*ncols+I-1]>voftol) || (C[J*ncols+I+1]<=voftol && C[J*ncols+I-1]<=voftol))
	{
	//Horz B below
			sr=C[J*ncols+I];
			sl=C[J*ncols+I];
			st=0.0;
			sb=1.0;
			if(v[J*ncols+I]>0.0)
			{
				if(v[J*ncols+I]*dt<=sr*dy)
				{
				Ft[J*ncols+I]=0.0;
				}
				else
				{
				Ft[J*ncols+I]=(v[J*ncols+I]*dt-(1-sr)*dy)*dx;
				}
			}
			if(v[(J-1)*ncols+I]<0.0)
			{			
				if(abs(v[(J-1)*ncols+I])*dt>=sr*dy)
				{
				Fb[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
				else
				{
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx;
				}
			}
	}
	else if(C[J*ncols+I+1]<=voftol && C[J*ncols+I-1]>voftol)
	{
	//Corn BL-B
		st=0.0;
		sr=0.0;
		sb=sqrt(C[J*ncols+I]*C[(J-1)*ncols+I]/C[J*ncols+I-1]);
		sl=C[J*ncols+I]/sb;
		if(sb>0.0 && sb<1.0 && sl>1.0)
		{
			sl=1.0;
			sb=C[J*ncols+I];
			st=C[J*ncols+I];
		}
		else if(sb>1.0)
		{
			sb=1.0;
			sl=C[J*ncols+I];
			sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt>(1-sl)*dy)
			{
			Ft[J*ncols+I]=(u[J*ncols+I]*dt-(1.0-sl)*dy)*sb*dx;
			}
			else
			{
			Ft[J*ncols+I]=0.0;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<sl*dy)
			{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*sb*dx;
			}
			else
			{
			Fb[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
		
	}
	else if(C[J*ncols+I+1]>voftol && C[J*ncols+I-1]<=voftol)
	{
	//Corn BR-B
		st=0.0;
		sl=0.0;
		sb=sqrt(C[J*ncols+I]*C[(J-1)*ncols+I]/C[J*ncols+I+1]);
		sr=C[J*ncols+I]/sb;
		if(sb>0.0 && sb<1.0 && sr>1.0)
		{
                    sr=1.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(sb>1.0)
		{
                    sb=1.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt>(1.0-sr)*dy)
			{
			Ft[J*ncols+I]=(v[J*ncols+I]*dt-(1.0-sr)*dy)*sb*dx;
			}
			else
			{
			Ft[J*ncols+I]=0.0;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<sr*dy)
			{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*sb*dx;
			}
			else
			{
			Fb[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
	}
	}
	else if(C[(J-1)*ncols+I]<=voftol  && C[(J+1)*ncols+I]>voftol)
	{
	if((C[J*ncols+I+1]>voftol && C[J*ncols+I-1]>voftol) || (C[J*ncols+I+1]<=voftol && C[J*ncols+I-1]<=voftol))
	{
	//Horz B above
		sr=C[J*ncols+I];
		sl=C[J*ncols+I];
		st=1.0;
		sb=0.0;
			if(v[J*ncols+I]>0.0)
			{
				if(v[J*ncols+I]*dt<sr*dy)
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
				if(abs(v[(J-1)*ncols+I])*dt>(1-sr)*dy)
				{
					Fb[J*ncols+I]=(abs(v[(J-1)*ncols+I])*dt-(1-sr)*dy)*dx;
				}
				else
				{
					Fb[J*ncols+I]=0.0;
				}
			}
	}
	else if(C[J*ncols+I+1]<=voftol && C[J*ncols+I-1]>voftol)
	{
		//Corn TL-B
		sb=0.0;
		sr=0.0;
		st=sqrt(C[J*ncols+I]*C[(J+1)*ncols+I]/C[J*ncols+I-1]);
		sl=C[J*ncols+I]/st;
		if(st>0.0 && st<1.0 && sl>1.0)
		{
                    sl=1.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(st>1.0)
		{
                    st=1.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<sl*dy)
			{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
			}
			else
			{
			Ft[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt>(1.0-sl)*dy)
			{
			Fb[J*ncols+I]=(abs(v[(J-1)*ncols+I])*dt-(1.0-sl)*dy)*st*dx;
			}
			else
			{
			Fb[J*ncols+I]=0.0;
			}
		}
	}
	else if(C[J*ncols+I+1]>voftol && C[J*ncols+I-1]<=voftol)
	{
	//Corn TR-B
		sb=0.0;
		sl=0.0;
		st=sqrt(C[J*ncols+I]*C[(J+1)*ncols+I]/C[J*ncols+I+1]);
		sr=C[J*ncols+I]/sb;
		if(st>0.0 && st<1.0 && sr>1.0)
		{
                    sr=1.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(st>1.0)
		{
                    st=1.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<sr*dy)
			{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
			}
			else
			{
			Ft[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
		if(v[J*ncols+I-1]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt>(1-sr)*dy)
			{
			Fb[J*ncols+I]=(abs(v[(J-1)*ncols+I])*dt-(1.0-sr)*dy)*st*dx;
			}
			else
			{
			Fb[J*ncols+I]=0.0;
			}
		}
	}
	}
	else if(C[(J-1)*ncols+I-1]>=(1.0-voftol)  && C[(J+1)*ncols+I+1]<(1.0-voftol))
	{
	if((C[J*ncols+I+1]>=(1.0-voftol) && C[J*ncols+I-1]>=(1.0-voftol)) || (C[J*ncols+I+1]<(1.0-voftol) && C[J*ncols+I-1]<(1.0-voftol)))
	{
	//Horz B below
			sr=C[J*ncols+I];
			sl=C[J*ncols+I];
			st=0.0;
			sb=1.0;
			if(v[J*ncols+I]>0.0)
			{
				if(v[J*ncols+I]*dt>(1-sr)*dy)
				{
					Ft[J*ncols+I]=(v[J*ncols+I]*dt-(1-sr)*dy)*dx;
				}
				else
				{
					Ft[J*ncols+I]=0.0;
				}
			}
			if(v[(J-1)*ncols+I]<0.0)
			{	
				if(abs(v[(J-1)*ncols+I])*dt<sr*dy)
				{
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx;
				}
				else
				{
					Fb[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
	
	}
	else if(C[J*ncols+I+1]>=(1.0-voftol) && C[J*ncols+I-1]<(1.0-voftol))
	{
	//Corn TL-W
		sb=1.0;
		sr=1.0;
		st=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J+1)*ncols+I]))/C[J*ncols+I-1]))))/2;
		sl=1.0-((1.0-C[J*ncols+I])/(1.0-st));
		if(st>0.0 && st<1.0 && sl<0.0)
		{
                    sl=0.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(st<0.0)
		{
                    st=0.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<(1.0-sl)*dy)
			{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
			}
			else
			{
			Ft[J*ncols+I]=st*dx*(1.0-sl)*dy+(v[J*ncols+I]*dt-(1.0-sl)*dy)*dx;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<sl*dx)
			{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx;
			}
			else
			{
			Fb[J*ncols+I]=sl*dx*dy+(abs(v[(J-1)*ncols+I])*dt-sl*dy)*st*dx;
			}
		}

	}
	else if(C[(J+1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]>=(1.0-voftol))
	{
	//Corn TR-W
		sb=1.0;
		sl=1.0;
		st=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J+1)*ncols+I]))/C[J*ncols+I+1]))))/2;
		sr=1.0-((1.0-C[J*ncols+I])/(1.0-st));
		if(st>0.0 && st<1.0 && sr<0.0)
		{
                    sr=0.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(st<0.0)
		{
                    st=0.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<(1.0-sr)*dy)
			{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
			}
			else
			{
			Ft[J*ncols+I]=st*dx*(1.0-sr)*dy+(v[J*ncols+I]*dt-(1.0-sr)*dy)*dx;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<sr*dy)
			{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*dx;
			}
			else
			{
			Fb[J*ncols+I]=sr*dy*dx+(abs(v[(J-1)*ncols+I])*dt-sr*dy)*st*dx;
			}
		}

	}
	}
	else if(C[J*ncols+I-1]<(1.0-voftol)  && C[J*ncols+I+1]>=(1.0-voftol))
	{
	if((C[(J+1)*ncols+I]>=(1.0-voftol) && C[(J-1)*ncols+I]>=(1.0-voftol)) || (C[(J+1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]<(1.0-voftol)))
	{
	//Vert B on right
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		sr=1.0;
		sl=0.0;
			if(v[J*ncols+I]>0.0)
			{
				Ft[J*ncols+I]=v[J*ncols+I]*dt*sb*dx;
			}
			if(v[(J-1)*ncols+I]<0.0)
			{			
				Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*sb*dx;
			}
	}
	else if(C[(J+1)*ncols+I]>=(1.0-voftol) && C[(J-1)*ncols+I]<(1.0-voftol))
	{
	//Corn BL-W
		st=1.0;
		sr=1.0;
		sb=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J-1)*ncols+I]))/C[J*ncols+I-1]))))/2;
		sl=1.0-((1.0-C[J*ncols+I])/(1.0-sb));
		if(sb>0.0 && sb<1.0 && sl<0.0)
		{
                    sl=0.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(sb<0.0)
		{
                    sb=0.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<sl*dy)
			{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*dx;
			}
			else
			{
			Ft[J*ncols+I]=sl*dy*dx+(v[J*ncols+I]*dt-sl*dy)*sb*dx;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<(1-sl)*dy)
			{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*sb*dx;
			}
			else
			{
			Fb[J*ncols+I]=sb*dx*(1.0-sl)*dy+(abs(v[(J-1)*ncols+I])*dt-(1.0-sl)*dy)*dx;
			}
		}

	}
	else if(C[(J+1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]>=(1.0-voftol))
	{
	//Corn BR-W
		st=1.0;
		sl=1.0;
		sb=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J-1)*ncols+I]))/C[J*ncols+I+1]))))/2;
		sr=1.0-((1.0-C[J*ncols+I])/(1.0-sb));
		if(sb>0.0 && sb<1.0 && sr<0.0)
		{
                    sr=0.0;
                    sb=C[J*ncols+I];
                    st=C[J*ncols+I];
		}
        else if(sb<0.0)
		{
                    sb=0.0;
                    sl=C[J*ncols+I];
                    sr=C[J*ncols+I];
		}
		if(v[J*ncols+I]>0.0)
		{
			if(v[J*ncols+I]*dt<sr*dy)
			{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*dx;
			}
			else
			{
			Ft[J*ncols+I]=sr*dx*dy+(v[J*ncols+I]*dt-sr*dy)*sb*dx;
			}
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			if(abs(v[(J-1)*ncols+I])*dt<(1.0-sr)*dy)
			{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*sb*dx;
			}
			else
			{
			Fb[J*ncols+I]=sb*dx*(1.0-sr)*dy+(abs(v[(J-1)*ncols+I])*dt-(1.0-sr)*dy)*dx;
			}
		}	
	}
	}
	
	else if(C[J*ncols+I-1]<=voftol  && C[J*ncols+I+1]<=voftol)
	{
	//Black finger
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		sr=0.0;
		sl=0.0;
		if(v[J*ncols+I]>0.0)
		{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*st*dx;
		}
		
	}
	else if(C[J*ncols+I-1]>=(1.0-voftol)  && C[J*ncols+I+1]>=(1.0-voftol))
	{
	//White finger
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		sr=1.0;
		sl=1.0;
		if(v[J*ncols+I]>0.0)
		{
			Ft[J*ncols+I]=v[J*ncols+I]*dt*st*dx;
		}
		if(v[(J-1)*ncols+I]<0.0)
		{
			Fb[J*ncols+I]=abs(v[(J-1)*ncols+I])*dt*st*dx;
		}	
	}
}
