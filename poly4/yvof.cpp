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
void FreeSurfaceCalc::RunYVOF()
{
float nx,ny,beta,alpha,st,sb,sr,sl,temp,temp2;
if(initialdirectionbool==0)
{
	//XSWEEP then YSWEEP
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	if((C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1])==0.0)
	{
		//Vertical Interface
		if(C[j*ncols+i-1]>C[j*ncols+i+1])
		{ 	//Vert B on left
			st=C[j*ncols+i];
			sb=C[j*ncols+i];
			sr=0.0;
			sl=1.0;
			if(u[j*ncols+i]>0.0)
			{
				if(u[j*ncols+i]*dt>(1.0-sb)*dx)
				{
				Fr[j*ncols+i]=(u[j*ncols+i]*dt-(1.0-sb)*dx)*dy;
				}
				else
				{
				Fr[j*ncols+i]=0.0;
				}
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				if(abs(u[j*ncols+i-1])*dt<sb*dx)
				{
				Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*dy;
				}
				else
				{
				Fl[j*ncols+i]=C[j*ncols+i]*dx*dy;
				}
			}
		}
		else if(C[j*ncols+i-1]<C[j*ncols+i+1])
		{
			//Vert B on right
		st=C[j*ncols+i];
		sb=C[j*ncols+i];
		sr=1.0;
		sl=0.0;
			if(u[j*ncols+i]>0.0)
			{
				if(u[j*ncols+i]*dt<sb*dx)
				{
				Fr[j*ncols+i]=u[j*ncols+i]*dt*dy;
				}
				else
				{
				Fr[j*ncols+i]=C[j*ncols+i]*dx*dy;
				}
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				if(abs(u[j*ncols+i-1])*dt>(1.0-sb)*dx)
				{
				Fl[j*ncols+i]=(abs(u[j*ncols+i-1])*dt-(1-sb)*dx)*dy;
				}
				else
				{
				Fl[j*ncols+i]=0.0;
				}
			}
		}
		else
		{
			cout<<"Error:1000 Vertical interface, Volume Fractions either side are equal"<<endl;
		}
	}
	else if((C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1])==0.0)
	{
		//Horizontal Interface
		if(C[(j-1)*ncols+i]>C[(j+1)*ncols+i])
		{ 	//Horiz B below
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=0.0;
			sb=1.0;
			if(u[j*ncols+i]>0.0)
			{
				Fr[j*ncols+i]=u[j*ncols+i]*dt*sr*dy;
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*sl*dy;
			}
		}
		else if(C[(j-1)*ncols+i]<C[(j+1)*ncols+i])
		{
			//Horz B above
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=1.0;
			sb=0.0;
			if(u[j*ncols+i]>0.0)
			{
				Fr[j*ncols+i]=u[j*ncols+i]*dt*sr*dy;
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*sl*dy;
			}
		}
		else
		{
			cout<<"Error:1000 with Horizontal Interface"<<endl;
		}
	}
	else
	{
	nx=(1.0/dx)*(C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1]);
	ny=(1.0/dy)*(C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1]);
	beta=atan(-nx/ny);
	alpha=atan((dx/dy)*tan(beta));
	if(alpha>0 && alpha<(PI/2))
	{
		//Either A or C
		if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//A
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			FluxCalcXsweepA(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//C
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=sb;
			temp2=sl;
			sb=st;
			sl=sr;
			st=temp;
			sr=temp2;
			FluxCalcXsweepC(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:1000 with 0<alpha<pi/2, A OR C"<<endl;
		}
	}
	else if(alpha>PI/2 && alpha<PI)
	{
    alpha=alpha-PI/2;
    beta=PI-beta;
        if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcXsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcXsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error: 1000 with pi/2<alpha<pi, D OR B"<<endl;
		}
	}
	else if(alpha<0 && alpha>(-PI/2))
	{
	alpha=(PI/2)-abs(alpha);
    beta=abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcXsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcXsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else
		{
			cout<<"Error:1000 with -pi/2<alpha<0, D OR B"<<endl;
		}
	}
	else if(alpha<(-PI/2) && alpha>(-PI))
	{
	alpha=PI-abs(alpha);
    beta=PI-abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//A (C')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		FluxCalcXsweepA(i,j,alpha,beta,st,sr,sb,sl);
	}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//C (A')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		temp=st;
        st=sb;
        sb=temp;
        temp2=sr;
        sr=sl;
        sl=temp2;
		FluxCalcXsweepC(i,j,alpha,beta,st,sr,sb,sl);
	}
	else
		{
			cout<<"Error:1000 with -pi<alpha<-pi/2, A OR C"<<endl;
		}
	}
	else
		{
			cout<<"Error:1000 Angle not within given ranges"<<endl; 
		}
	}
	}
	else if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true)
	{
		//Full Cells
		if(u[j*ncols+i]>0.0)
		{
		Fr[j*ncols+i]=u[j*ncols+i]*dt*dy;
		}
		/*else
		{
			cout<<"Eror:1000 No outward velocity out right face"<<endl;
		}*/
		if(u[j*ncols+i-1]<0.0)
		{
		Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*dy;
		}
		/*else
		{
			cout<<"Error:1000 No outward velocity out the left face"<<endl;
		}*/
	}
	/*else
	{
		cout<<"Error:1000 Not a CV cell or doesnt have any fluid in it"<<endl; 
	}*/
	}
}
ImposePeriodicFluxBC();
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true)
	{
	C[j*ncols+i]=(C[j*ncols+i]+((Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy)))/(1.0-(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]));
	if(C[j*ncols+i]>1.0)
	{
		//cout<<"Error:1000, C is greater than 1"<<endl;
		if(u[j*ncols+i]>0.0 && u[j*ncols+i-1]>0.0)
		{
		Fr[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]<0.0)
		{
		Fl[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]>0.0)
		{
			int n=0;
			if(v[j*ncols+i]>0.0 && v[j*ncols+i-1]>0.0)
			{
			Fr[(j+1)*ncols+i-1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(v[j*ncols+i]<0.0 && v[j*ncols+i-1]<0.0)
			{
			Fr[(j-1)*ncols+i-1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(v[j*ncols+i]<0.0 && v[j*ncols+i-1]>0.0)
			{
				cout<<"Error:1000, C is greater than 1, all velocities point inward"<<endl;
				C[j*ncols+i]=1.0;
			}
		}
	}
	if(C[j*ncols+i]<0.0)
	{
		//cout<<"Error:1000, C is less than 0"<<endl;
		C[j*ncols+i]=0.0;
	}
	}
	/*else
	{
		cout<<"Error:1000 Not a CV cell"<<endl;
	}*/
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	if((C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1])==0.0)
	{
		if(C[j*ncols+i-1]>C[j*ncols+i+1])
		{
		//Vert B on left
			st=C[j*ncols+i];
			sb=C[j*ncols+i];
			sr=0.0;
			sl=1.0;
			if(v[j*ncols+i]>0.0)
			{
				Ft[j*ncols+i]=v[j*ncols+i]*dt*sb*dx;
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*sb*dx;
			}
		}
		else if(C[j*ncols+i-1]<C[j*ncols+i+1])
		{
			//Vert B on right
			st=C[j*ncols+i];
			sb=C[j*ncols+i];
			sr=1.0;
			sl=0.0;
			if(v[j*ncols+i]>0.0)
			{
				Ft[j*ncols+i]=v[j*ncols+i]*dt*sb*dx;
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*sb*dx;
			}
		}
		else
		{
			cout<<"Error:0100 Vertical Interface"<<endl;
		}
	}
	else if((C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1])==0.0)
	{
		if(C[(j-1)*ncols+i]>C[(j+1)*ncols+i])
		{
		//Horz B below
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=0.0;
			sb=1.0;
			if(v[j*ncols+i]>0.0)
			{
				if(v[j*ncols+i]*dt>(1-sr)*dy)
				{
					Ft[j*ncols+i]=(v[j*ncols+i]*dt-(1-sr)*dy)*dx;
				}
				else
				{
					Ft[j*ncols+i]=0;
				}
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				if(abs(v[(j-1)*ncols+i])*dt<sr*dy)
				{
					Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*dx;//*****
				}
				else
				{
					Fb[j*ncols+i]=C[j*ncols+i]*dx*dy;
				}
			}
		}
		else if(C[(j-1)*ncols+i]<C[(j+1)*ncols+i])
		{
			//Horz B above
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=1.0;
			sb=0.0;
			if(v[j*ncols+i]>0.0)
			{
				if(v[j*ncols+i]*dt<sr*dy)
				{
					Ft[j*ncols+i]=v[j*ncols+i]*dt*dx;
				}
				else
				{
					Ft[j*ncols+i]=C[j*ncols+i]*dx*dy;;
				}
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				if(abs(v[(j-1)*ncols+i])*dt>(1-sr)*dy)
				{
					Fb[j*ncols+i]=(abs(v[(j-1)*ncols+i])*dt-(1-sr)*dy)*dx;//******
				}
				else
				{
					Fb[j*ncols+i]=0.0;
				}
			}
		}
		else
		{
			cout<<"Error:0100 Horizontal interface"<<endl;
		}
	}
	else
	{
	nx=(1.0/dx)*(C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1]);
	ny=(1.0/dy)*(C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1]);
	beta=atan(-nx/ny);
	alpha=atan((dx/dy)*tan(beta));
	if(alpha>0 && alpha<(PI/2))
	{
		//Either A or C
		if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//A
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			FluxCalcYsweepA(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//C
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=sb;
			temp2=sl;
			sb=st;
			sl=sr;
			st=temp;
			sr=temp2;
			FluxCalcYsweepC(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:0100 alpha>0 && alpha<(PI/2), A or C"<<endl;
		}
	}
	else if(alpha>PI/2 && alpha<PI)
	{
    alpha=alpha-PI/2;
    beta=PI-beta;
        if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcYsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcYsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:0100 alpha>PI/2 && alpha<PI, D or B"<<endl;
		}
	}
	else if(alpha<0 && alpha>(-PI/2))
	{
	alpha=(PI/2)-abs(alpha);
    beta=abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcYsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcYsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else
	{
		cout<<"Error:0100 alpha<0 && alpha>(-PI/2), D or B"<<endl;
	}
	}
	else if(alpha<(-PI/2) && alpha>(-PI))
	{
	alpha=PI-abs(alpha);
    beta=PI-abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//A (C')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		FluxCalcYsweepA(i,j,alpha,beta,st,sr,sb,sl);
	}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//C (A')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		temp=st;
        st=sb;
        sb=temp;
        temp2=sr;
        sr=sl;
        sl=temp2;
		FluxCalcYsweepC(i,j,alpha,beta,st,sr,sb,sl);
	}
	else
	{
		cout<<"Error:0100 alpha<(-PI/2) && alpha>(-PI), A or C"<<endl;
	}
	}
	else
	{
		cout<<"Error:0100 Alpha not within range"<<endl;
	}
	}
	}
	else if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true)
	{
		if(v[j*ncols+i]>0.0)
		{
		Ft[j*ncols+i]=v[j*ncols+i]*dt*dx;
		}
		/*else
		{
			cout<<"Error: 0100 no outward velocity at top boundary"<<endl;
		}*/
		if(v[(j-1)*ncols+i]<0.0)
		{
		Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*dx;
		}
		/*else
		{
			cout<<"Error:0100 outward velocity at bottom bounadry"<<endl;
		}*/
	}
	/*else
	{
		cout<<"Error:0100 Not A CV cell"<<endl;
	}*/
	}
}
ImposePeriodicFluxBC();
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true)
	{
	//C[j*ncols+i]=C[j*ncols+i]*(1.0+(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]))+(Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy);
	C[j*ncols+i]=(C[j*ncols+i]+((Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy)))/(1.0-(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]));
	if(C[j*ncols+i]>1.0)
	{
		//cout<<"Error:0100, C is greater than 1"<<endl;
		if(v[j*ncols+i]>0.0 && v[(j-1)*ncols+i]>0.0)
		{
		Ft[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(v[j*ncols+i]<0.0 && v[(j-1)*ncols+i]<0.0)
		{
		Fb[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(v[j*ncols+i]<0.0 && v[(j-1)*ncols+i]>0.0)
		{
			int n=0;
			if(u[j*ncols+i]>0.0 && u[j*ncols+i-1]>0.0)
			{
			Ft[(j-1)*ncols+i+1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]<0.0)
			{
			Ft[(j-1)*ncols+i-1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]>0.0)
			{
				cout<<"Error:0100, C is greater than 1, all velocities point inward"<<endl;
				C[j*ncols+i]=1.0;
			}
		}
	}
	if(C[j*ncols+i]<0.0)
	{
		//cout<<"Error:0100, C is less than 0"<<endl;
		C[j*ncols+i]=0.0;
	}
	}
	/*else
	{
		cout<<"Error:0100 Not A CV cell"<<endl;
	}*/
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		Fr[j*ncols+i]=0.0;
		Fl[j*ncols+i]=0.0;
		Ft[j*ncols+i]=0.0;
		Fb[j*ncols+i]=0.0;
	}
}
initialdirectionbool=1;
}
else if(initialdirectionbool==1)
{
	//YSWEEP then XSWEEP
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	if((C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1])==0.0)
	{
		if(C[j*ncols+i-1]>C[j*ncols+i+1])
		{
		//Vert B on left
			st=C[j*ncols+i];
			sb=C[j*ncols+i];
			sr=0.0;
			sl=1.0;
			if(v[j*ncols+i]>0.0)
			{
				Ft[j*ncols+i]=v[j*ncols+i]*dt*sb*dx;
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*sb*dx;
			}
		}
		else if(C[j*ncols+i-1]<C[j*ncols+i+1])
		{
			//Vert B on right
			st=C[j*ncols+i];
			sb=C[j*ncols+i];
			sr=1.0;
			sl=0.0;
			if(v[j*ncols+i]>0.0)
			{
				Ft[j*ncols+i]=v[j*ncols+i]*dt*sb*dx;
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*sb*dx;
			}
		}
		else
		{
			cout<<"Error:0010 vertical boundary"<<endl;
		}
	}
	else if((C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1])==0.0)
	{
		if(C[(j-1)*ncols+i]>C[(j+1)*ncols+i])
		{
		//Horz B below
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=0.0;
			sb=1.0;
			if(v[j*ncols+i]>0.0)
			{
				if(v[j*ncols+i]*dt>(1-sr)*dy)
				{
					Ft[j*ncols+i]=(v[j*ncols+i]*dt-(1-sr)*dy)*dx;
				}
				else
				{
					Ft[j*ncols+i]=0;
				}
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				if(abs(v[(j-1)*ncols+i])*dt<sr*dy)
				{
					Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*dx;//*****
				}
				else
				{
					Fb[j*ncols+i]=C[j*ncols+i]*dx*dy;
				}
			}
			}
		else if(C[(j-1)*ncols+i]<C[(j+1)*ncols+i])
		{
			//Horz B above
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=1.0;
			sb=0.0;
			if(v[j*ncols+i]>0.0)
			{
				if(v[j*ncols+i]*dt<sr*dy)
				{
					Ft[j*ncols+i]=v[j*ncols+i]*dt*dx;
				}
				else
				{
					Ft[j*ncols+i]=C[j*ncols+i]*dx*dy;;
				}
			}
			if(v[(j-1)*ncols+i]<0.0)
			{			
				if(abs(v[(j-1)*ncols+i])*dt>(1-sr)*dy)
				{
					Fb[j*ncols+i]=(abs(v[(j-1)*ncols+i])*dt-(1-sr)*dy)*dx;//*****
				}
				else
				{
					Fb[j*ncols+i]=0.0;
				}
			}
		}
		else
		{
			cout<<"Error:0010 Horizontal interface"<<endl;
		}
	}
	else
	{
	nx=(1.0/dx)*(C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1]);
	ny=(1.0/dy)*(C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1]);
	beta=atan(-nx/ny);
	alpha=atan((dx/dy)*tan(beta));
	if(alpha>0 && alpha<(PI/2))
	{
		//Either A or C
		if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//A
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			FluxCalcYsweepA(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//C
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=sb;
			temp2=sl;
			sb=st;
			sl=sr;
			st=temp;
			sr=temp2;
			FluxCalcYsweepC(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:0010 alpha>0 && alpha<(PI/2), A or C"<<endl;
		}
	}
	else if(alpha>PI/2 && alpha<PI)
	{
    alpha=alpha-PI/2;
    beta=PI-beta;
        if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcYsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcYsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:0010 alpha>PI/2 && alpha<PI, B or D"<<endl;
		}
	}
	else if(alpha<0 && alpha>(-PI/2))
	{
	alpha=(PI/2)-abs(alpha);
    beta=abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcYsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcYsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else
	{
		cout<<"Error:0010 alpha<0 && alpha>(-PI/2), D or B"<<endl;
	}
	}
	else if(alpha<(-PI/2) && alpha>(-PI))
	{
	alpha=PI-abs(alpha);
    beta=PI-abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//A (C')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		FluxCalcYsweepA(i,j,alpha,beta,st,sr,sb,sl);
	}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//C (A')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		temp=st;
        st=sb;
        sb=temp;
        temp2=sr;
        sr=sl;
        sl=temp2;
		FluxCalcYsweepC(i,j,alpha,beta,st,sr,sb,sl);
	}
	else
	{
		cout<<"Error:0010 alpha<(-PI/2) && alpha>(-PI), A or C"<<endl;
	}
	}
	else
	{
		cout<<"Error:0010 Alpha out of range"<<endl;
	}
	}
	}
	else if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true)
	{
		if(v[j*ncols+i]>0.0)
		{
		Ft[j*ncols+i]=v[j*ncols+i]*dt*dx;
		}
		/*else
		{
			cout<<"Error:0010 No outward velocity out of the top boundary"<<endl;
		}*/
		if(v[(j-1)*ncols+i]<0.0)
		{
		Fb[j*ncols+i]=abs(v[(j-1)*ncols+i])*dt*dx;
		}
		/*else
		{
			cout<<"Error:0010 No outward velocity out of the bottom boundary"<<endl;
		}*/
	}
	/*else
	{
		cout<<"Error:0010 Not a CV cell"<<endl;
	}*/
	}
}
ImposePeriodicFluxBC();
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true)
	{
	C[j*ncols+i]=(C[j*ncols+i]+((Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy)))/(1.0-(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]));
	if(C[j*ncols+i]>1.0)
	{
		cout<<"Error:0010, C is greater than 1"<<endl;
				if(v[j*ncols+i]>0.0 && v[(j-1)*ncols+i]>0.0)
		{
		Ft[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(v[j*ncols+i]<0.0 && v[(j-1)*ncols+i]<0.0)
		{
		Fb[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(v[j*ncols+i]<0.0 && v[(j-1)*ncols+i]>0.0)
		{
			int n=0;
			if(u[j*ncols+i]>0.0 && u[j*ncols+i-1]>0.0)
			{
			Ft[(j-1)*ncols+i+1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]<0.0)
			{
			Ft[(j-1)*ncols+i-1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]>0.0)
			{
				cout<<"Error:0010, C is greater than 1, all velocities point inward"<<endl;
				C[j*ncols+i]=1.0;
			}
		}
	}
	if(C[j*ncols+i]<0.0)
	{
		cout<<"Error:0010, C is less than 0"<<endl;
		C[j*ncols+i]=0.0;
	}
	}
	/*else
	{
		cout<<"Error:0010 Not a CV cell"<<endl;
	}*/
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	if((C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1])==0.0)
	{
		//Vertical Interface
		if(C[j*ncols+i-1]>C[j*ncols+i+1])
		{ 	//Vert B on left
			st=C[j*ncols+i];
			sb=C[j*ncols+i];
			sr=0.0;
			sl=1.0;
			if(u[j*ncols+i]>0.0)
			{
				if(u[j*ncols+i]*dt>(1.0-sb)*dx)
				{
				Fr[j*ncols+i]=(u[j*ncols+i]*dt-(1.0-sb)*dx)*dy;
				}
				else
				{
				Fr[j*ncols+i]=0.0;
				}
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				if(abs(u[j*ncols+i-1])*dt<sb*dx)
				{
				Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*dy;
				}
				else
				{
				Fl[j*ncols+i]=C[j*ncols+i]*dx*dy;
				}
			}
		}
		else if(C[j*ncols+i-1]<C[j*ncols+i+1])
		{
			//Vert B on right
		st=C[j*ncols+i];
		sb=C[j*ncols+i];
		sr=1.0;
		sl=0.0;
			if(u[j*ncols+i]>0.0)
			{
				if(u[j*ncols+i]*dt<sb*dx)
				{
				Fr[j*ncols+i]=u[j*ncols+i]*dt*dy;
				}
				else
				{
				Fr[j*ncols+i]=C[j*ncols+i]*dx*dy;
				}
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				if(abs(u[j*ncols+i-1])*dt>(1.0-sb)*dx)
				{
				Fl[j*ncols+i]=(abs(u[j*ncols+i-1])*dt-(1-sb)*dx)*dy;
				}
				else
				{
				Fl[j*ncols+i]=0.0;
				}
			}
		}
		else
		{
			cout<<"Error:0001 Vertical interface"<<endl;
		}
	}
	else if((C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1])==0.0)
	{
		//Horizontal Interface
		if(C[(j-1)*ncols+i]>C[(j+1)*ncols+i])
		{ 	//Horiz B below
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=0.0;
			sb=1.0;
			if(u[j*ncols+i]>0.0)
			{
				Fr[j*ncols+i]=u[j*ncols+i]*dt*sr*dy;
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*sl*dy;
			}
		}
		else if(C[(j-1)*ncols+i]<C[(j+1)*ncols+i])
		{
			//Horz B above
			sr=C[j*ncols+i];
			sl=C[j*ncols+i];
			st=1.0;
			sb=0.0;
			if(u[j*ncols+i]>0.0)
			{
				Fr[j*ncols+i]=u[j*ncols+i]*dt*sr*dy;
			}
			if(u[j*ncols+i-1]<0.0)
			{			
				Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*sl*dy;
			}
		}
		else
		{
			cout<<"Error:0001 Horizontal Interface"<<endl;
		}
	}
	else
	{
	nx=(1.0/dx)*(C[(j+1)*ncols+i+1]+2.0*C[j*ncols+i+1]+C[(j-1)*ncols+i+1]-C[(j+1)*ncols+i-1]-2.0*C[j*ncols+i-1]-C[(j-1)*ncols+i-1]);
	ny=(1.0/dy)*(C[(j+1)*ncols+i+1]+2.0*C[(j+1)*ncols+i]+C[(j+1)*ncols+i-1]-C[(j-1)*ncols+i+1]-2.0*C[(j-1)*ncols+i]-C[(j-1)*ncols+i-1]);
	beta=atan(-nx/ny);
	alpha=atan((dx/dy)*tan(beta));
	if(alpha>0 && alpha<(PI/2))
	{
		//Either A or C
		if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//A
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			FluxCalcXsweepA(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//C
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=sb;
			temp2=sl;
			sb=st;
			sl=sr;
			st=temp;
			sr=temp2;
			FluxCalcXsweepC(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:0001 alpha>0 && alpha<(PI/2), A or C"<<endl;
		}
	}
	else if(alpha>PI/2 && alpha<PI)
	{
    alpha=alpha-PI/2;
    beta=PI-beta;
        if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcXsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcXsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
		else
		{
			cout<<"Error:0001 alpha>PI/2 && alpha<PI, D or B"<<endl;
		}
	}
	else if(alpha<0 && alpha>(-PI/2))
	{
	alpha=(PI/2)-abs(alpha);
    beta=abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//D (C')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sr;
			sr=sb;
			sb=sl;
			sl=temp;
			//Switch II and III
			FluxCalcXsweepDIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
		{
			//B (A')
			YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
			temp=st;
			st=sl;
			sl=sb;
			sb=sr;
			sr=temp;
			//Switch II and III
			FluxCalcXsweepBIIrevIII(i,j,alpha,beta,st,sr,sb,sl);
		}
	else
	{
		cout<<"Error:0001 alpha<0 && alpha>(-PI/2), D or B"<<endl;
	}
	}
	else if(alpha<(-PI/2) && alpha>(-PI))
	{
	alpha=PI-abs(alpha);
    beta=PI-abs(beta);
	if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])>(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])<(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//A (C')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		FluxCalcXsweepA(i,j,alpha,beta,st,sr,sb,sl);
	}
	else if ((C[(j+1)*ncols+i+1]+2*C[j*ncols+i+1]+C[(j-1)*ncols+i+1])<(C[(j+1)*ncols+i-1]+2*C[j*ncols+i-1]+C[(j-1)*ncols+i-1])||(C[(j+1)*ncols+i-1]+2*C[(j+1)*ncols+i]+C[(j+1)*ncols+i+1])>(C[(j-1)*ncols+i-1]+2*C[(j-1)*ncols+i]+C[(j-1)*ncols+i+1]))
	{
		//C (A')
		YVOFInterfaceConstantCalc(i,j,alpha,st,sr,sb,sl);
		temp=st;
        st=sb;
        sb=temp;
        temp2=sr;
        sr=sl;
        sl=temp2;
		FluxCalcXsweepC(i,j,alpha,beta,st,sr,sb,sl);
	}
	else
	{
		cout<<"Error:0001 alpha<(-PI/2) && alpha>(-PI), A or C"<<endl;
	}
	}
	else
	{
		cout<<"Error:0001 alpha not within range"<<endl;
	}
	}
	}
	else if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true)
	{
		//Full Cells
		if(u[j*ncols+i]>0.0)
		{
		Fr[j*ncols+i]=u[j*ncols+i]*dt*dy;
		}
		/*else
		{
			cout<<"Error:0001 No outward velocity out of the right hand cell boundary"<<endl;
		}*/
		if(u[j*ncols+i-1]<0.0)
		{
		Fl[j*ncols+i]=abs(u[j*ncols+i-1])*dt*dy;
		}
		/*else
		{
			cout<<"Error:0001 No outward velocity out of the left hand cell boundary"<<endl;
		}*/
	}
	/*else
	{
		cout<<"Error:0001 Cell is not in the CV"<<endl;
	}*/
	}
}
ImposePeriodicFluxBC();
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true)
	{
	//C[j*ncols+i]=C[j*ncols+i]*(1.0+(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]))+(Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy);
	C[j*ncols+i]=(C[j*ncols+i]+((Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy)))/(1.0-(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]));
	if(C[j*ncols+i]>1.0)
	{
		//cout<<"Error:0001, C is greater than 1"<<endl;
		if(u[j*ncols+i]>0.0 && u[j*ncols+i-1]>0.0)
		{
		Fr[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]<0.0)
		{
		Fl[j*ncols+i]+=(C[j*ncols+i]-1.0)*dx*dy;
		C[j*ncols+i]=1.0;
		}
		else if(u[j*ncols+i]<0.0 && u[j*ncols+i-1]>0.0)
		{
			int n=0;
			if(v[j*ncols+i]>0.0 && v[j*ncols+i-1]>0.0)
			{
			Fr[(j+1)*ncols+i-1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(v[j*ncols+i]<0.0 && v[j*ncols+i-1]<0.0)
			{
			Fr[(j-1)*ncols+i-1]+=(C[j*ncols+i]-1.0)*dx*dy;
			C[j*ncols+i]=1.0;
			}
			else if(v[j*ncols+i]<0.0 && v[j*ncols+i-1]>0.0)
			{
				cout<<"Error:0001, C is greater than 1, all velocities point inward"<<endl;
				C[j*ncols+i]=1.0;
			}
		}
	}
	if(C[j*ncols+i]<0.0)
	{
		//cout<<"Error:0001, C is less than 0"<<endl;
		C[j*ncols+i]=0.0;
	}
	}
	/*else
	{
		cout<<"Error:0001 Cell is not in the CV"<<endl;
	}*/
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		Fr[j*ncols+i]=0.0;
		Fl[j*ncols+i]=0.0;
		Ft[j*ncols+i]=0.0;
		Fb[j*ncols+i]=0.0;
	}
}
initialdirectionbool=0;
}
}

void FreeSurfaceCalc::YVOFInterfaceConstantCalc(int &I, int &J, float &ALPHA, float &ST, float &SR, float &SB, float &SL)
{
if(ALPHA<(PI/4))
{
	if(C[J*ncols+I]<=(0.5*tan(ALPHA)))
	{
	//case I
	ST=0.0;
	SR=sqrt(2.0*C[J*ncols+I]*tan(ALPHA));
	SB=sqrt(2.0*C[J*ncols+I]*(1.0/tan(ALPHA)));
	SL=0.0;
	}
	else if(C[J*ncols+I]<=(1.0-0.5*tan(ALPHA)))
	{
	//Case II
	ST=0.0;
	SR=C[J*ncols+I]+0.5*tan(ALPHA);
	SB=1.0;
	SL=C[J*ncols+I]-0.5*tan(ALPHA);
	}
	else
	{
	//Case IV
	ST=1.0-sqrt(2.0*(1.0-C[J*ncols+I])*(1.0/tan(ALPHA)));
	SR=1.0;
	SB=1.0;
	SL=1.0-sqrt(2.0*(1.0-C[J*ncols+I])*tan(ALPHA));
	}
}
else
{
	if(C[J*ncols+I]<=(0.5*(1.0/tan(ALPHA))))
	{
	//case I
	ST=0.0;
	SR=sqrt(2.0*C[J*ncols+I]*tan(ALPHA));
	SB=sqrt(2.0*C[J*ncols+I]*(1.0/tan(ALPHA)));
	SL=0.0;
	}
	else if(C[J*ncols+I]<=(1.0-0.5*(1.0/tan(ALPHA))))
	{
	//Case III
	ST=C[J*ncols+I]-0.5*(1.0/tan(ALPHA));
	SR=1.0;
	SB=C[J*ncols+I]+0.5*(1.0/tan(ALPHA));
	SL=1.0;
	}
	else
	{
	//Case IV
	ST=1.0-sqrt(2.0*(1.0-C[J*ncols+I])*(1.0/tan(ALPHA)));
	SR=1.0;
	SB=1.0;
	SL=1.0-sqrt(2.0*(1.0-C[J*ncols+I])*tan(ALPHA));
	}
}
}
void FreeSurfaceCalc::ImposePeriodicFluxBC()
{
	//Boundary Cells on far right equal boundary cells on far left, ie u[j*ncols+0]=u[j*ncols+ncols-2], with free slip at the vertical boundaries. 
	//Lower and upper Boundary Cells are still considered rigid, no slip.
	for(int j=0;j<nrows;j++)
	{
	Fr[j*ncols+0]=Fr[j*ncols+ncols-2];
	Fl[j*ncols+ncols-1]=Fl[j*ncols+1];
	v[j*ncols+0]=v[j*ncols+ncols-2];
	v[j*ncols+ncols-1]=v[j*ncols+1];
	}
	for(int i=0;i<ncols;i++)
	{
	Ft[0*ncols+i]=0.0;
	Ft[(nrows-2)*ncols+i]=0.0;
	Fb[(nrows-1)*ncols+i]=0.0;
	Fb[1*ncols+i]=0.0;
	}		
}