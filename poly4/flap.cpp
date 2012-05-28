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

void FreeSurfaceCalc::RunFLAP()
{
if(initialdirectionbool==0)
{
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
		{
		FLAPXsweep(i,j);
		}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	//C[j*ncols+i]=C[j*ncols+i]*(1.0+(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]))+(Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy);
	C[j*ncols+i]=(C[j*ncols+i]+(Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy))/(1.0-(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]));
	}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
		{
		FLAPYsweep(i,j);
		}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	C[j*ncols+i]=C[j*ncols+i]*(1.0+(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]))+(Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy);
	//C[j*ncols+i]=(C[j*ncols+i]+(Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy))/(1.0-(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]));
	}
	}
}
initialdirectionbool=1;
}
else
{
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
		{
		FLAPYsweep(i,j);
		}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	//C[j*ncols+i]=C[j*ncols+i]*(1.0+(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]))+(Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy);
	C[j*ncols+i]=(C[j*ncols+i]+(Fb[(j+1)*ncols+i]+Ft[(j-1)*ncols+i]-Fb[j*ncols+i]-Ft[j*ncols+i])/(dx*dy))/(1.0-(dt/dy)*(v[j*ncols+i]-v[(j-1)*ncols+i]));
	}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
		if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
		{
		FLAPXsweep(i,j);
		}
	}
}
for(int j=0;j<nrows;j++)
{
	for(int i=0;i<ncols;i++)
	{
	if(fl[j*ncols+i].C==true && fl[j*ncols+i].F==true && fl[j*ncols+i].S==true)
	{
	//C[j*ncols+i]=(C[j*ncols+i]+(Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy))/(1.0-(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]));
	C[j*ncols+i]=C[j*ncols+i]*(1.0+(dt/dx)*(u[j*ncols+i]-u[j*ncols+i-1]))+(Fl[j*ncols+i+1]+Fr[j*ncols+i-1]-Fl[j*ncols+i]-Fr[j*ncols+i])/(dx*dy);
	}
	}
}
initialdirectionbool=0;
}
}

