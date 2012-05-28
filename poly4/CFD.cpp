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
#define g 9.81
using namespace std;
//Constructor
CFD::CFD(int inncols, int innrows, float indx, float indy, float inRe, float inrho, float indynvis, bool inSMACbool, bool inSLICbool, bool inFLAPbool, bool inYVOFbool):ncols(inncols),nrows(innrows), dx(indx), dy(indy), Re(inRe), rho(inrho), dynvis(indynvis), SMACbool(inSMACbool), SLICbool(inSLICbool), FLAPbool(inFLAPbool), YVOFbool(inYVOFbool)
{
}
CFD::~CFD()
{
}
//Methods

void CFD::SetInitialConditions()
{
string str;
str="InitialU.dat";
SetVecofFloats(str,u);
str="InitialV.dat";
SetVecofFloats(str,v);
str="InitialP.dat";
SetVecofFloats(str,p);

if(SMACbool==1)
{
str="InitialX.dat";
SetVecofFloats(str,x);
str="InitialY.dat";
SetVecofFloats(str,y);
vector<int> mask;
str="mask.dat";
SetVecofInts(str,mask);
int maskncols=ncols+2;
int masknrows=nrows+2;
Convertflag(mask,maskncols,masknrows,fl);
}
else if(SLICbool==1 || FLAPbool==1 || YVOFbool==1)
{
vector<int> mask;
str="mask.dat";
SetVecofInts(str,mask);
int maskncols=ncols+2;
int masknrows=nrows+2;
Convertflag(mask,maskncols,masknrows,fl);
//str="InitialC.dat";
//SetVecofFloats(str,C);
SetInitialVelocityandPressure();
FreeSurfaceConditions();
FindInterfaceandSetC();
SetVolumeFractionInBoundary();
int minitial=0;
float tinitial=0.0;
PrintData(minitial,tinitial);
}
}
void CFD::SetVecofFloats(string myStr,vector<float> &myVec)
{
ifstream file(myStr);
while (file)
{
	string line;
	getline(file, line);
	istringstream is(line);
	while (!is.eof())
	{
	float data;
	is >> data;
	myVec.push_back(data);
	}
	myVec.pop_back();
}
}
void CFD::SetVecofInts(string myStr,vector<int> &myVec)
{
ifstream file(myStr);
while (file)
{
	string line;
	getline(file, line);
	istringstream is(line);
	while (!is.eof())
	{
	int data;
	is >> data;
	myVec.push_back(data);
	}
	myVec.pop_back();
}
}
void CFD::Convertflag(vector<int> &maskvec, int &ncol, int &nrow, flaglist &FL)
{
for(int j=1;j<=nrow-2;j++){
	for(int i=1;i<=ncol-2;i++){
		if(maskvec[j*ncol+i]==1)
		{
		//flag f1(0,0,0);
		//FL.push_back(f1);
			if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==1)
			{
				//CV borders below
				flag f1(0,0,0,1,0,0,0);
				FL.push_back(f1);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==1)
			{
				//CV borders above
				flag f2(0,0,0,0,1,0,0);
				FL.push_back(f2);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==0)
			{
				//CV borders to the right
				flag f3(0,1,0,0,0,0,0);
				FL.push_back(f3);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==1)
			{
				//CV borders to the left
				flag f4(0,0,1,0,0,0,0);
				FL.push_back(f4);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==0)
			{
				//CV borders below and to the right
				flag f5(0,1,0,1,0,0,0);
				FL.push_back(f5);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==1)
			{
				//CV borders below and to the left
				flag f6(0,0,1,1,0,0,0);
				FL.push_back(f6);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==0)
			{
				//CV borders above and to the right
				flag f7(0,1,0,0,1,0,0);
				FL.push_back(f7);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==1)
			{
				//CV borders above and to the left
				flag f8(0,0,1,0,1,0,0);
				FL.push_back(f8);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==1)
			{
				//Boundary cell that doesnt border CV
				flag f9(0,0,0,0,0,0,0);
				FL.push_back(f9);
			}
			else
			{
				cout<<"Non Viable Boundary cell orientation at i:"<<i<<"j: "<<j<<endl;
			}
		}
		else
		{
		//flag f2(1,0,0);
		//FL.push_back(f2);
			if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==0)
			{
				//Boundary borders below
				flag fl1(1,1,1,0,1,0,0);
				FL.push_back(fl1);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==0)
			{
				//Boundary borders above
				flag fl2(1,1,1,1,0,0,0);
				FL.push_back(fl2);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==1)
			{
				//Boundary borders to the right
				flag fl3(1,0,1,1,1,0,0);
				FL.push_back(fl3);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==0)
			{
				//Boundary borders to the left
				flag fl4(1,1,0,1,1,0,0);
				FL.push_back(fl4);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==1)
			{
				//Boundary borders below and to the right
				flag fl5(1,0,1,0,1,0,0);
				FL.push_back(fl5);
			}
			else if(maskvec[(j-1)*ncol+i]==1 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==0)
			{
				//Boundary borders below and to the left
				flag fl6(1,1,0,0,1,0,0);
				FL.push_back(fl6);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==0 && maskvec[j*ncol+i+1]==1)
			{
				//Boundary borders above and to the right
				flag fl7(1,0,1,1,0,0,0);
				FL.push_back(fl7);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==1 && maskvec[j*ncol+i-1]==1 && maskvec[j*ncol+i+1]==0)
			{
				//Boundary borders above and to the left
				flag fl8(1,1,0,1,0,0,0);
				FL.push_back(fl8);
			}
			else if(maskvec[(j-1)*ncol+i]==0 && maskvec[(j+1)*ncol+i]==0 && maskvec[j*ncol+i-1]==0.0 && maskvec[j*ncol+i+1]==0)
			{
				//CV that doesnt border Boundary cell
				flag fl9(1,1,1,1,1,0,0);
				FL.push_back(fl9);
			}
			else
			{
				cout<<"ERROR: Non Viable CV cell orientation at i: "<<i<<" j: "<<j<<". Exit"<<endl;
				//return;
			}
		}
	}
}

}
void CFD::PrintData(int &myMval, float &myTval)
{
{
char bufferu[50];
int filenameu=sprintf(bufferu, "Ut000%d.dat", (myMval));
char bufferv[50];
int filenamev=sprintf(bufferv, "Vt000%d.dat", (myMval));
char bufferp[50];
int filenamep=sprintf(bufferp, "Pt000%d.dat", (myMval));
char bufferf[50];
int filenamef=sprintf(bufferf, "Ft000%d.dat", (myMval));
char bufferC[50];
int filenamec=sprintf(bufferC, "Ct000%d.dat", (myMval));
ofstream foutu(bufferu,ios::binary);
ofstream foutv(bufferv,ios::binary);
ofstream foutp(bufferp,ios::binary);
ofstream foutfluidflag(bufferf,ios::binary);
ofstream foutC(bufferC,ios::binary);
//ofstream foutomega("Finalomega.dat",ios::binary);
//ofstream foutpsi("Finalpsi.dat",ios::binary);
for(int j=0;j<=(nrows-1);j++)
	{
		for(int i=0;i<=(ncols-1);i++)
		{
			foutu<<u[j*ncols+i]<<" ";
			foutv<<v[j*ncols+i]<<" ";
			foutp<<p[j*ncols+i]<<" ";
			foutfluidflag<<fl[j*ncols+i].F<<" ";
			foutC<<C[j*ncols+i]<<" ";
			//foutomega<<omega[j*ncols+i]<<",";
			//foutpsi<<psi[j*ncols+i]<<",";
		}
		foutu<<endl;
		foutv<<endl;
		foutp<<endl;
		foutfluidflag<<endl;
		foutC<<endl;
		//foutomega<<endl;
		//foutpsi<<endl;
	}
}
char bufferextra[50];
int filenamec=sprintf(bufferextra, "SuppInfo000%d.dat", (myMval));
ofstream foutextra(bufferextra,ios::binary);
foutextra<<"Re: "<<Re<<", m: "<<myMval<<", t: "<<myTval;//<<", L2normx: "<<L2normx<<", L2normy: "<<L2normy<<endl;
}
void CFD::VorticityCalc()
{

}
void CFD::StreamlineCalc()
{

}
void CFD::SetVolumeFractionInBoundary()
{
	for(int j=0;j<nrows;j++)
	{
		for(int i=0;i<ncols;i++)
		{
			if(fl[j*ncols+i].C==false)
			{
				C[j*ncols+i]=1.0;
			}
		}
	}
}

void CFD::SetInitialVelocityandPressure()
{
	float ytemp,jinterface;
	for(int j=1;j<=nrows-2;j++)
	{
		for(int i=1;i<=ncols-2;i++)
		{
			//Equations for a plane progressive wave
			//lambda=20m,d=18m,A=0.5 dt=PI/omega (or fractions there of), domain=2*lambda 
			int lambda=20,d=18;//meters
			int rho=1000;//kilograms per meter cubed
			float A=0.25;//meters
			int omega=10;//seconds
			int t=0;//Check wjether it makes a difference starting at t=0 or t=-pi/2omega
			float k=2*PI/lambda;
			ytemp=d+A*sin((2*PI/lambda)*i*dx-omega*t);
			jinterface=floor(ytemp/dy);
			if(j<jinterface)
			{
			u[j*ncols+i]=A*(2*PI/lambda)*sqrt((g*lambda)/(2*PI))*exp((2*PI/lambda)*(j*dy-d))*sin((2*PI/lambda)*i*dx-omega*t);
			v[j*ncols+i]=-A*(2*PI/lambda)*sqrt((g*lambda)/(2*PI))*exp((2*PI/lambda)*(j*dy-d))*cos((2*PI/lambda)*i*dx-omega*t);
			p[j*ncols+i]=rho*g*(A*sin((2*PI/lambda)*((i-1)*dx+0.5*dx)-omega*t)*exp((2*PI/lambda)*(j*dy-d))-j*dy-d);
			}
		}
	}
}
void CFD::FindInterfaceandSetC()
{
	float st,sb,sl,sr,ytemp,ytempplus,ytempminus,xtemp;
	float j,jplus,jminus;
	int lambda=20,d=18;//meters
	int rho=1000;//kilograms per meter cubed
	float A=0.25;//meters
	int omega=10;//seconds
	float k=2*PI/lambda;
	int t=0;//Check whether it makes a difference starting at t=0 or t=-pi/2omega
	float c1;
	for(int jiterator=0;jiterator<=nrows-1;jiterator++)
	{
		for(int i=0;i<=ncols-1;i++)
		{
			c1=0.0;
			C.push_back(c1);
		}
	}
	for(int i=0;i<=ncols-2;i++)//uncertain about upperlimit
	{
	ytemp=d+A*sin(k*i*dx-omega*t);
	j=floor(ytemp/dy);
	sl=fmod(ytemp,dy)/dy;
	ytempplus=d+A*sin(k*(i+1)*dx-omega*t);
	jplus=floor(ytempplus/dy);
	ytempminus=d+A*sin(k*(i-1)*dx-omega*t);
	jminus=floor(ytempminus/dy);
	//Note that jplus and jminus will not be required in a boundary cell
 	if(sl!=0.0)
	{
	if(fl[(j+1)*ncols+i+1].C==true && fl[(j+1)*ncols+i+1].S==false)
	{
	if(jplus!=j)
	{		
		if(jplus>j)
		{
			if((fmod(ytempplus,dy)/dy)==0)
			{
			st=0.0;
            sb=1.0;
            sr=1.0;
            C[(j+1)*ncols+i+1]=sl+(1-sl)/2.0;
			}
			else
			{
			xtemp=(asin(((j+1)*dy-d)/A)+omega*t)/k;
			while(xtemp<i*dx)
			{
				xtemp+=lambda/2.0;
			}
			sr=1.0;
			st=1.0-fmod(xtemp,dx)/dx;
			sb=1.0;
			C[(j+1)*ncols+i+1]=1.0-((1.0-sl)*(1.0-st)/2.0);
			}
		}
		else
		{
			xtemp=(asin((j*dy-d)/A)+omega*t)/k;
			while(xtemp<i*dx)
			{
				xtemp+=lambda/2.0;
			}
			sr=0.0;
			sb=1.0-fmod(xtemp,dx)/dx;
			st=0.0;
			C[(j+1)*ncols+i+1]=sl*sr/2.0;
		}
	}
	else
	{
		sr=fmod(ytempplus,dy)/dy;
		if(sr>sl)
		{
			C[(j+1)*ncols+i+1]=sl+((sr-sl)/2.0);
		}
		else
		{
			C[(j+1)*ncols+i+1]=sr+((sl-sr)/2.0);
		}
	}
	fl[(j+1)*ncols+i+1].S=true;
	}
	if(fl[(j+1)*ncols+i].C==true && fl[(j+1)*ncols+i].S==false)
	{
	sr=sl;
	if(jminus!=j)
	{
		if(jminus>j)
		{
			xtemp=(asin(((j+1)*dy-d)/A)+omega*t)/k;
			while(xtemp<(i-1)*dx)
			{
				xtemp+=lambda/2.0;
			}
			sl=1.0;
			st=fmod(xtemp,dx)/dx;
			sb=1.0;
			C[(j+1)*ncols+i]=1.0-((1.0-sr)*(1.0-st)/2.0);
		}
		else
		{
			xtemp=(asin((j*dy-d)/A)+omega*t)/k;
			while(xtemp<(i-1)*dx)
			{
				xtemp+=lambda/2.0;
			}
			sl=0.0;
			st=0.0;
			sb=1-fmod(xtemp,dx)/dx;
			C[(j+1)*ncols+i]=sb*sr/2.0;
		}
	}
	else
	{
		if(sl>sr)
		{
			C[(j+1)*ncols+i]=sr+((sl-sr)/2.0);
		}
		else
		{
			C[(j+1)*ncols+i]=sl+((sr-sl)/2);
		}
	}
	fl[(j+1)*ncols+i].S=true;
	}
	}
	else
	{

		//Interface goes through BL cell corner
		if(jplus!=j)
		{
			if(jplus>j)
			{
				if(fl[(j+1)*ncols+i+1].C==true && fl[(j+1)*ncols+i+1].S==false)
				{
				xtemp=(asin(((j+1)*dy-d)/A)+omega*t)/k;
				while(xtemp<i*dx)
				{
					xtemp+=lambda/2.0;
				}
				sr=1.0;
				st=fmod(xtemp,dx)/dx;
				sb=1.0;
				C[(j+1)*ncols+i+1]=st+((1.0-st)/2.0);
				fl[(j+1)*ncols+i+1].S=true;
				}
			}
			else if(jplus==j-1)
			{
				if(fl[j*ncols+i+1].C==true && fl[j*ncols+i+1].S==false)
				{
				sr=fmod(ytempplus,dy)/dy;
				sb=1.0;
				st=0.0;
				sl=1.0;
				C[j*ncols+i+1]=sr+((1.0-sr)/2.0);
				fl[j*ncols+i+1].S=true;
				}
			}
			else
			{
				if(fl[j*ncols+i+1].C==true && fl[j*ncols+i+1].S==false)
				{
				xtemp=(asin(((j-1)*dy-d)/A)+omega*t)/k;
				while(xtemp<i*dx)
				{
					xtemp+=lambda/2.0;
				}
				sb=fmod(xtemp,dx)/dx;
				sr=0.0;
				sl=1.0;
				st=0.0;
				C[j*ncols+i+1]=sb/2.0;
				fl[j*ncols+i+1].S=true;
				}
			}
		}
		else
		{
			if(fl[(j+1)*ncols+i+1].C==true && fl[(j+1)*ncols+i+1].S==false)
			{
			sb=1.0;
			st=0.0;
			sl=0.0;
			sr=fmod(ytempplus,dy)/dy;
			C[(j+1)*ncols+i+1]=sr/2.0;
			fl[(j+1)*ncols+i+1].S=true;
			}
		}
		if(jminus!=j)
		{
			if(jminus>j)
			{
			if(fl[(j+1)*ncols+i].C==true && fl[(j+1)*ncols+i].S==false)
			{
				xtemp=(asin(((j+1)*dy-d)/A)+omega*t)/k;
				while(xtemp<i*dx)
				{
					xtemp+=lambda/2.0;
				}
				st=fmod(xtemp,dx)/dx;
				sb=1.0;
				sl=1.0;
				sr=0.0;
				C[(j+1)*ncols+i]=st+((1.0-st)/2.0);
				fl[(j+1)*ncols+i].S=true;
			}
			}
			else if(jminus==j-1)
			{
							
			if(fl[j*ncols+i].C==true && fl[j*ncols+i].S==false)
			{
				sl=fmod(ytempminus,dy)/dy;
				sr=1.0;
				st=0.0;
				sb=1.0;
				C[j*ncols+i]=sl+((1.0-sl)/2.0);
				fl[j*ncols+i].S=true;
			}
			}
			else
			{
				if(fl[j*ncols+i].C==true && fl[j*ncols+i].S==false)
				{
				xtemp=(asin(((j-1)*dy-d)/A)+omega*t)/k;
				while(xtemp<i*dx)
				{
					xtemp+=lambda/2.0;
				}
				sb=1.0-fmod(xtemp,dx)/dx;
				sr=1.0;
				sl=0.0;
				st=0.0;
				C[j*ncols+i]=sb/2.0;
				fl[j*ncols+i].S=true;
				}
			}
		}
		else
		{
			if(fl[(j+1)*ncols+i].C==true && fl[(j+1)*ncols+i].S==false)
			{
			sl=fmod(ytempminus,dy)/dy;
			sr=0.0;
			st=0.0;
			sb=1.0;
			C[(j+1)*ncols+i]=sl/2.0;
			fl[(j+1)*ncols+i].S=true;
			}
		}
	}
}
int jmax;
for(int i=0;i<=ncols-2;i++)//uncertain about upperlimit
{
ytemp=d+A*sin(k*i*dx-omega*t);
jmax=floor(ytemp/dy);
	for(int j=1;j<=jmax;j++)
	{
		if(fl[j*ncols+i+1].C==true && fl[j*ncols+i+1].S==false)
		{
		C[j*ncols+i+1]=1.0;
		}
	}
}
}