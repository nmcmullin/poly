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
void CFD::FreeSurfaceConditions()
{
for(int j=0;j<nrows;j++)
{
		for(int i=0;i<ncols;i++)
		{
			if(fl[j*ncols+i].F==1 && fl[j*ncols+i].S==1)
			{
				if(fl[j*ncols+i].B==1 && fl[(j-1)*ncols+i].F==0)
				{
					if(fl[(j*ncols+i)].L==1 && fl[(j*ncols+i-1)].F==0)
					{
					p[j*ncols+i]=(1.0/(2.0*Re))*(((u[(j+1)*ncols+i]+u[(j+1)*ncols+i-1]-u[j*ncols+i]-u[j*ncols+i-1])/dy)+((v[j*ncols+i+1]+v[(j-1)*ncols+i+1]-v[j*ncols+i]-v[(j-1)*ncols+i])/dx));
					u[j*ncols+i-1]=u[j*ncols+i];
					v[(j-1)*ncols+i]=v[j*ncols+i];
					}
					else if(fl[(j*ncols+i)].R==1 && fl[(j*ncols+i+1)].F==0)
					{
					p[j*ncols+i]=(-1.0/(2.0*Re))*(((u[(j+1)*ncols+i]+u[(j+1)*ncols+i-1]-u[j*ncols+i]-u[j*ncols+i-1])/dy)+((v[j*ncols+i]+v[(j-1)*ncols+i]-v[j*ncols+i-1]-v[(j-1)*ncols+i-1])/dx));
					u[j*ncols+i]=u[j*ncols+i-1];
					v[(j-1)*ncols+i]=v[j*ncols+i];
					}
					else
					{
					p[j*ncols+i]=(2.0/Re)*((v[j*ncols+i]-v[(j-1)*ncols+i])/dy);
					u[(j-1)*ncols+i]=u[j*ncols+i]+(dy/dx)*(v[(j-1)*ncols+i+1]-v[(j-1)*ncols+i]);
					}
				}
				else if(fl[j*ncols+i].T==1 && fl[(j+1)*ncols+i].F==0)
				{
					if(fl[(j*ncols+i)].L==1 && fl[(j*ncols+i-1)].F==0)
					{
					p[j*ncols+i]=(-1.0/(2.0*Re))*(((u[j*ncols+i]+u[j*ncols+i-1]-u[(j-1)*ncols+i]-u[(j-1)*ncols+i-1])/dy)+((v[j*ncols+i+1]+v[(j-1)*ncols+i+1]-v[j*ncols+i]-v[(j-1)*ncols+i])/dx));
					u[j*ncols+i-1]=u[j*ncols+i];
					v[(j-1)*ncols+i]=v[j*ncols+i];
					}
					else if(fl[(j*ncols+i)].R==1 && fl[(j*ncols+i+1)].F==0)
					{
					p[j*ncols+i]=(1.0/(2.0*Re))*(((u[j*ncols+i]+u[j*ncols+i-1]-u[(j-1)*ncols+i]-u[(j-1)*ncols+i-1])/dy)+((v[j*ncols+i]+v[(j-1)*ncols+i]-v[j*ncols+i-1]-v[(j-1)*ncols+i-1])/dx));
					u[j*ncols+i]=u[j*ncols+i-1];
					v[j*ncols+i]=v[(j-1)*ncols+i];
					}
					else
					{
					p[j*ncols+i]=(2.0/Re)*((v[j*ncols+i]-v[(j-1)*ncols+i])/dy);
					u[(j+1)*ncols+i-1]=u[j*ncols+i-1]-(dy/dx)*(v[j*ncols+i]-v[j*ncols+i-1]);
					}
				}
				else if(fl[j*ncols+i].L==1 && fl[j*ncols+i-1].F==0)
				{
				p[j*ncols+i]=(2.0/Re)*((u[j*ncols+i]-u[j*ncols+i-1])/dx);
				v[(j-1)*ncols+i-1]=v[(j-1)*ncols+i]+(dx/dy)*(u[j*ncols+i-1]-u[(j-1)*ncols+i-1]);
				}
				else if(fl[j*ncols+i].R==1 && fl[j*ncols+i+1].F==0)
				{
				p[j*ncols+i]=(2.0/Re)*((u[j*ncols+i]-u[j*ncols+i-1])/dx);
				v[j*ncols+i+1]=v[j*ncols+i]-(dx/dy)*(u[(j+1)*ncols+i]-u[j*ncols+i]);
				}
			}
		}
}

/*			if(fl[j*ncols+i+1].F==1 && fl[j*ncols+i-1].F==1 && fl[(j+1)*ncols+i].F==1 && fl[(j-1)*ncols+i].F==0 && fl[j*ncols+i].B==1)
			{
				//Fluid above, Absence of fluid below
				p[j*ncols+i]=(2.0/Re)*((v[j*ncols+i]-v[(j-1)*ncols+i])/dy);
				u[(j-1)*ncols+i]=u[j*ncols+i]+(dy/dx)*(v[(j-1)*ncols+i+1]-v[(j-1)*ncols+i]);
			}
			else if(fl[j*ncols+i+1].F==1 && fl[j*ncols+i-1].F==1 && fl[(j+1)*ncols+i].F==0 && fl[(j-1)*ncols+i].F==1 && fl[j*ncols+i].T==1)
			{
				//Fluid beneath, Absence of fluid above
				p[j*ncols+i]=(2.0/Re)*((v[j*ncols+i]-v[(j-1)*ncols+i])/dy);
				u[(j+1)*ncols+i-1]=u[j*ncols+i-1]-(dy/dx)*(v[j*ncols+i]-v[j*ncols+i-1]);
			}
			else if(fl[j*ncols+i+1].F==1 && fl[j*ncols+i-1].F==0 && fl[(j+1)*ncols+i].F==1 && fl[(j-1)*ncols+i].F==1 && fl[j*ncols+i].L==1)
			{
				//Fluid to the right, Absence of fluid left
				p[j*ncols+i]=(2.0/Re)*((u[j*ncols+i]-u[j*ncols+i-1])/dx);
				v[(j-1)*ncols+i-1]=v[(j-1)*ncols+i]+(dx/dy)*(u[j*ncols+i-1]-u[(j-1)*ncols+i-1]);
			}
			else if(fl[j*ncols+i+1].F==0 && fl[j*ncols+i-1].F==1 && fl[(j+1)*ncols+i].F==1 && fl[(j-1)*ncols+i].F==1 && fl[j*ncols+i].R==1)
			{
				//Fluid to the left, Absence of fluid right
				p[j*ncols+i]=(2.0/Re)*((u[j*ncols+i]-u[j*ncols+i-1])/dx);
				v[j*ncols+i+1]=v[j*ncols+i]-(dx/dy)*(u[(j+1)*ncols+i]-u[j*ncols+i]);
			}
			else if(fl[j*ncols+i+1].F==0 && fl[j*ncols+i-1].F==1 && fl[(j+1)*ncols+i].F==0 && fl[(j-1)*ncols+i].F==1 && fl[j*ncols+i].T==1 && fl[j*ncols+i].R==1)
			{
				//Fluid in bottom left, Absence of fluid above and to the right
				p[j*ncols+i]=(1.0/(2.0*Re))*(((u[j*ncols+i]+u[j*ncols+i-1]-u[(j-1)*ncols+i]-u[(j-1)*ncols+i-1])/dy)+((v[j*ncols+i]+v[(j-1)*ncols+i]-v[j*ncols+i-1]-v[(j-1)*ncols+i-1])/dx));
				u[j*ncols+i]=u[j*ncols+i-1];
				v[j*ncols+i]=v[(j-1)*ncols+i];
			}
			else if(fl[j*ncols+i+1].F==1 && fl[j*ncols+i-1].F==0 && fl[(j+1)*ncols+i].F==0 && fl[(j-1)*ncols+i].F==1 && fl[j*ncols+i].T==1 && fl[j*ncols+i].L==1)
			{
				//Fluid in bottom right, Absence of fluid above and to the left
				p[j*ncols+i]=(-1.0/(2.0*Re))*(((u[j*ncols+i]+u[j*ncols+i-1]-u[(j-1)*ncols+i]-u[(j-1)*ncols+i-1])/dy)+((v[j*ncols+i+1]+v[(j-1)*ncols+i+1]-v[j*ncols+i]-v[(j-1)*ncols+i])/dx));
				u[j*ncols+i-1]=u[j*ncols+i];
				v[(j-1)*ncols+i]=v[j*ncols+i];
			}
			else if(fl[j*ncols+i+1].F==0 && fl[j*ncols+i-1].F==1 && fl[(j+1)*ncols+i].F==1 && fl[(j-1)*ncols+i].F==0 && fl[j*ncols+i].B==1 && fl[j*ncols+i].R==1)
			{
				//Fluid in top left, Absence of fluid below and to the right
				p[j*ncols+i]=(-1.0/(2.0*Re))*(((u[(j+1)*ncols+i]+u[(j+1)*ncols+i-1]-u[j*ncols+i]-u[j*ncols+i-1])/dy)+((v[j*ncols+i]+v[(j-1)*ncols+i]-v[j*ncols+i-1]-v[(j-1)*ncols+i-1])/dx));
				u[j*ncols+i]=u[j*ncols+i-1];
				v[(j-1)*ncols+i]=v[j*ncols+i];
			}
			else if(fl[j*ncols+i+1].F==1 && fl[j*ncols+i-1].F==0 && fl[(j+1)*ncols+i].F==1 && fl[(j-1)*ncols+i].F==0 && fl[j*ncols+i].B==1 && fl[j*ncols+i].L==1)
			{
				//Fluid in top right, Absence of fluid below and to the left
				p[j*ncols+i]=(1.0/(2.0*Re))*(((u[(j+1)*ncols+i]+u[(j+1)*ncols+i-1]-u[j*ncols+i]-u[j*ncols+i-1])/dy)+((v[j*ncols+i+1]+v[(j-1)*ncols+i+1]-v[j*ncols+i]-v[(j-1)*ncols+i])/dx));
				u[j*ncols+i-1]=u[j*ncols+i];
				v[(j-1)*ncols+i]=v[j*ncols+i];
			}
			}
		}
}*/
}