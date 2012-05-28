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

void FreeSurfaceCalc::FLAPXsweep(int &I, int &J)
{
	float sr,sl,st,sb;
	if((C[J*ncols+I-1]<(1.0-voftol) && C[J*ncols+I-1]>voftol) && (C[J*ncols+I+1]<(1.0-voftol) && C[J*ncols+I+1]>voftol))//(1.0-C[j*ncols+i])<(1.0-voftol))//i.e if the cell contains an interface
	{
		//Type I, horizontal interface with the interface at a height of the VF*dy of the lower fluid (the one with VF=1 below) 
		sr=C[J*ncols+I];
		sl=C[J*ncols+I];
		if(C[(J-1)*ncols+I]>=C[(J+1)*ncols+I])
		{
			st=0.0;
			sb=1.0;
			if(u[J*ncols+I]>0.0)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*sr*dy;
			}
			if(u[J*ncols+I-1]<0.0)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*sr*dy;
			}
		}
		else
		{
			st=1.0;
			sb=0.0;
			if(u[J*ncols+I]>0.0)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*sr*dy;
			}
			if(u[J*ncols+I-1]<0.0)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*sr*dy;
			}	
		}
	}
	else if(C[J*ncols+I-1]>voftol  && C[J*ncols+I+1]<=voftol)
	{
	if((C[(J+1)*ncols+I]>voftol && C[(J-1)*ncols+I]>voftol) || (C[(J+1)*ncols+I]<=voftol && C[(J-1)*ncols+I]<=voftol))
	{
	//Vert B on left
			st=C[J*ncols+I];
			sb=C[J*ncols+I];
			sr=0.0;
			sl=1.0;
			if(u[J*ncols+I]>0.0)
			{
				if(u[J*ncols+I]*dt>(1.0-sb)*dx)
				{
				Fr[J*ncols+I]=(u[J*ncols+I]*dt-(1.0-sb)*dx)*dy;
				}
				else
				{
				Fr[J*ncols+I]=0.0;
				}
			}
			if(u[J*ncols+I-1]<0.0)
			{			
				if(abs(u[J*ncols+I-1])*dt<sb*dx)
				{
				Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
				}
				else
				{
				Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
	}
	else if(C[(J+1)*ncols+I]<=voftol && C[(J-1)*ncols+I]>voftol)
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt>(1-sb)*dx)
			{
			Fr[J*ncols+I]=(u[J*ncols+I]*dt-(1.0-sb)*dx)*sl*dy;
			}
			else
			{
			Fr[J*ncols+I]=0.0;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<sb*dx)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*sl*dy;
			}
			else
			{
			Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
		
	}
	else if(C[(J+1)*ncols+I]>voftol && C[(J-1)*ncols+I]<=voftol)
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt>(1-st)*dx)
			{
			Fr[J*ncols+I]=(u[J*ncols+I]*dt-(1.0-st)*dx)*sl*dy;
			}
			else
			{
			Fr[J*ncols+I]=0.0;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<st*dx)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*sl*dy;
			}
			else
			{
			Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}	
	}
	}
	else if(C[J*ncols+I-1]<=voftol  && C[J*ncols+I+1]>voftol)
	{
	if((C[(J+1)*ncols+I]>voftol && C[(J-1)*ncols+I]>voftol) || (C[(J+1)*ncols+I]<=voftol && C[(J-1)*ncols+I]<=voftol))
	{
	//Vert B on Right
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		sr=1.0;
		sl=0.0;
			if(u[J*ncols+I]>0.0)
			{
				if(u[J*ncols+I]*dt<sb*dx)
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
				if(abs(u[J*ncols+I-1])*dt>(1.0-sb)*dx)
				{
				Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-(1-sb)*dx)*dy;
				}
				else
				{
				Fl[J*ncols+I]=0.0;
				}
			}
	}
	else if(C[(J+1)*ncols+I]<=voftol && C[(J-1)*ncols+I]>voftol)
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<sb*dx)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*sr*dy;
			}
			else
			{
			Fr[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt>(1-sb)*dx)
			{
			Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-(1.0-sb)*dx)*sr*dy;
			}
			else
			{
			Fl[J*ncols+I]=0.0;
			}
		}
	}
	else if(C[(J+1)*ncols+I]>voftol && C[(J-1)*ncols+I]<=voftol)
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<st*dx)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*sr*dy;
			}
			else
			{
			Fr[J*ncols+I]=C[J*ncols+I]*dx*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt>(1-st)*dx)
			{
			Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-(1.0-st)*dx)*sr*dy;
			}
			else
			{
			Fl[J*ncols+I]=0.0;
			}
		}
	}
	}
	else if(C[J*ncols+I-1]>=(1.0-voftol)  && C[J*ncols+I+1]<(1.0-voftol))
	{
	if((C[(J+1)*ncols+I]>=(1.0-voftol) && C[(J-1)*ncols+I]>=(1.0-voftol)) || (C[(J+1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]<(1.0-voftol)))
	{
	//Vert B on left
			st=C[J*ncols+I];
			sb=C[J*ncols+I];
			sr=0.0;
			sl=1.0;
			if(u[J*ncols+I]>0.0)
			{
				if(u[J*ncols+I]*dt>(1.0-sb)*dx)
				{
				Fr[J*ncols+I]=(u[J*ncols+I]*dt-(1.0-sb)*dx)*dy;
				}
				else
				{
				Fr[J*ncols+I]=0.0;
				}
			}
			if(u[J*ncols+I-1]<0.0)
			{			
				if(abs(u[J*ncols+I-1])*dt<sb*dx)
				{
				Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
				}
				else
				{
				Fl[J*ncols+I]=C[J*ncols+I]*dx*dy;
				}
			}
	
	}
	else if(C[(J+1)*ncols+I]>=(1.0-voftol) && C[(J-1)*ncols+I]<(1.0-voftol))
	{
	//Corn BR-W
		st=1.0;
		sl=1.0;
		sb=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J-1)*ncols+I]))/(1.0-C[J*ncols+I+1])))))/2;
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<(1.0-sb)*dx)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*sr*dy;
			}
			else
			{
			Fr[J*ncols+I]=sr*dy*(1.0-sb)*dx+(u[J*ncols+I]*dt-(1.0-sb)*dx)*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<sb*dx)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
			}
			else
			{
			Fl[J*ncols+I]=sb*dx*dy+(abs(u[J*ncols+I-1])*dt-sb*dx)*sr*dy;
			}
		}

	}
	else if(C[(J+1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]>=(1.0-voftol))
	{
	//Corn TR-W
		sb=1.0;
		sl=1.0;
		st=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J+1)*ncols+I]))/(1.0-C[J*ncols+I+1])))))/2;
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<(1.0-st)*dx)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*sr*dy;
			}
			else
			{
			Fr[J*ncols+I]=sr*dy*(1.0-st)*dx+(u[J*ncols+I]*dt-(1.0-st)*dx)*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<st*dx)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
			}
			else
			{
			Fl[J*ncols+I]=st*dx*dy+(abs(u[J*ncols+I-1])*dt-st*dx)*sr*dy;
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
			if(u[J*ncols+I]>0.0)
			{
				if(u[J*ncols+I]*dt<sb*dx)
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
				if(abs(u[J*ncols+I-1])*dt>(1.0-sb)*dx)
				{
				Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-(1-sb)*dx)*dy;
				}
				else
				{
				Fl[J*ncols+I]=0.0;
				}
			}
	}
	else if(C[(J+1)*ncols+I]>=(1.0-voftol) && C[(J-1)*ncols+I]<(1.0-voftol))
	{
	//Corn BL-W
		st=1.0;
		sr=1.0;
		sb=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J-1)*ncols+I]))/(1.0-C[J*ncols+I-1])))))/2;
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<sb*dx)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*dy;
			}
			else
			{
			Fr[J*ncols+I]=sb*dx*dy+(u[J*ncols+I]*dt-sb*dx)*sl*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<(1-sb)*dx)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*sl*dy;
			}
			else
			{
			Fl[J*ncols+I]=sl*dy*(1.0-sb)*dx+(abs(u[J*ncols+I-1])*dt-(1.0-sb)*dx)*dy;
			}
		}

	}
	else if(C[(J+1)*ncols+I]<(1.0-voftol) && C[(J-1)*ncols+I]>=(1.0-voftol))
	{
	//Corn TL-W
		sb=1.0;
		sr=1.0;
		st=1.0-(sqrt(4.0-4.0*(1.0-(((1.0-C[J*ncols+I])*(1.0-C[(J+1)*ncols+I]))/(1.0-C[J*ncols+I-1])))))/2;
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
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<st*dx)
			{
			Fr[J*ncols+I]=u[J*ncols+I]*dt*dy;
			}
			else
			{
			Fr[J*ncols+I]=st*dx*dy+(u[J*ncols+I]*dt-st*dx)*sl*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<(1-st)*dx)
			{
			Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*sl*dy;
			}
			else
			{
			Fl[J*ncols+I]=sl*dy*(1.0-st)*dx+(abs(u[J*ncols+I-1])*dt-(1.0-st)*dx)*dy;
			}
		}
	}
	}
	else if(C[J*ncols+I-1]<=voftol  && C[J*ncols+I+1]<=voftol)
	{
	//Vert Black finger
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		sr=0.0;
		sl=0.0;
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<0.5*(1.0-st)*dx)
			{
				Fr[J*ncols+I]=0.0;
			}
			else if(u[J*ncols+I]*dt<(0.5*(1.0-st)*dx+st*dx))
			{
				Fr[J*ncols+I]=(u[J*ncols+I]*dt-0.5*(1.0-st)*dx)*dy;
			}
			else
			{
				Fr[J*ncols+I]=st*dx*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<0.5*(1.0-st)*dx)
			{
				Fl[J*ncols+I]=0.0;
			}
			else if(abs(u[J*ncols+I-1])*dt<(0.5*(1.0-st)*dx+st*dx))
			{
				Fl[J*ncols+I]=(abs(u[J*ncols+I-1])*dt-0.5*(1.0-st)*dx)*dy;
			}
			else
			{
				Fl[J*ncols+I]=st*dx*dy;
			}
		}
		
	}
	else if(C[J*ncols+I-1]>=(1.0-voftol)  && C[J*ncols+I+1]>=(1.0-voftol))
	{
	//Vert White finger
		st=C[J*ncols+I];
		sb=C[J*ncols+I];
		sr=1.0;
		sl=1.0;
		if(u[J*ncols+I]>0.0)
		{
			if(u[J*ncols+I]*dt<0.5*st*dx)
			{
				Fr[J*ncols+I]=u[J*ncols+I]*dt*dy;
			}
			else if(u[J*ncols+I]*dt<(0.5*st*dx+(1.0-st)*dx))
			{
				Fr[J*ncols+I]=0.5*st*dx*dy;
			}
			else
			{
				Fr[J*ncols+I]=0.5*st*dx*dy+(u[J*ncols+I]*dt-(0.5*st*dx+(1.0-st)*dx))*dy;
			}
		}
		if(u[J*ncols+I-1]<0.0)
		{
			if(abs(u[J*ncols+I-1])*dt<0.5*st*dx)
			{
				Fl[J*ncols+I]=abs(u[J*ncols+I-1])*dt*dy;
			}
			else if(abs(u[J*ncols+I-1])*dt<(0.5*st*dx+(1.0-st)*dx))
			{
				Fl[J*ncols+I]=0.5*st*dx*dy;
			}
			else
			{
				Fl[J*ncols+I]=0.5*st*dx*dy+(abs(u[J*ncols+I-1])*dt-(0.5*st*dx+(1.0-st)*dx))*dy;
			}
		}	
	}
}