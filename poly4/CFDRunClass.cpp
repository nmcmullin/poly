#include "CFDRunClass.h"
#include "NavierStokesCalc.h"
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
//constructor
CFDRunClass::CFDRunClass(int inncols,int innrows,float indx, float indy,float inRe, float inrho, float indynvis, float intmax, float inmaxNumIt, float inveleps, float inpeps, float ingx, float ingy, float inw, float intsf, float ingsf, bool inGSbool, bool inCGbool, bool inSMACbool, int innpart, bool inSLICbool, bool inFLAPbool, bool inYVOFbool, float invoftol): CFD(inncols,innrows,indx,indy,inRe,inrho,indynvis,inSMACbool,inSLICbool,inFLAPbool,inYVOFbool),NavierStokesCalc(inncols,innrows,indx,indy,inRe,inrho,indynvis,inSMACbool,inSLICbool,inFLAPbool,inYVOFbool,intmax,inmaxNumIt,inveleps,inpeps,ingx,ingy,inw,intsf,ingsf,inGSbool,inCGbool),FreeSurfaceCalc(inncols,innrows,indx,indy,inRe,inrho,indynvis,inSMACbool,innpart,inSLICbool,inFLAPbool,inYVOFbool,invoftol)
{
}
CFDRunClass::~CFDRunClass()
{
}
//Methods
void CFDRunClass::Run()
{
SetInitialConditions();
bool velcheck;
UpdateFluidFlag();
SetPreliminaryGuessValues();
float t=0;
int m=0;
do
{
 	RunNavierStokes();
	velcheck=UpdateVelocity();
	RunFreeSurfaceCalc();
	t+=dt;
	m+=1;
	vector<float> timevec;
	timevec.push_back(0.02);timevec.push_back(0.05);timevec.push_back(0.1);timevec.push_back(0.2);timevec.push_back(0.3);timevec.push_back(0.4);timevec.push_back(0.5);timevec.push_back(0.6);timevec.push_back(0.7);timevec.push_back(0.8);timevec.push_back(0.9);timevec.push_back(1.0);
	//timevec.push_back(0.01);timevec.push_back(0.02);timevec.push_back(0.05);timevec.push_back(0.1);timevec.push_back(0.2);timevec.push_back(0.3);timevec.push_back(0.4);
	for(int mit=0;mit<=11;mit++)
	{
	if(t==timevec[mit] || ((t-dt)<timevec[mit] && t>timevec[mit]))
	{
	PrintData(m,t);
	PrintFreeSurfaceData(m);
	}
	}
}while((t<tmax)&&(velcheck==true));
PrintData(m,t);
}