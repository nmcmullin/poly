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
using namespace std;
#ifndef CFD_H
#define CFD_H
class CFD
{
public:
	CFD(int inncols,int innrows, float indx, float indy, float inRe, float inrho, float indynvis, bool inSMACbool,bool inSLICbool,bool inFLAPbool,bool inYVOFbool);
	virtual ~CFD();
	virtual void SetInitialConditions();
	virtual void VorticityCalc();
	virtual void StreamlineCalc();
	virtual void PrintData(int &myMval, float &myTval);
	void FreeSurfaceConditions();
protected:
	int const ncols,nrows;
	float dt;
	float const dx,dy,Re,rho,dynvis;
	vector<float> u,v,p,x,y,C;
	flaglist fl;
	bool SMACbool,SLICbool,FLAPbool,YVOFbool;
private:
	void SetVecofFloats(string myStr,vector<float> &myVec);
	void SetVecofInts(string myStr,vector<int> &myVec);
	void Convertflag(vector<int> &maskvec, int &ncol, int &nrow, flaglist &FL);
	void SetVolumeFractionInBoundary();
	void FindInterfaceandSetC();
	void SetInitialVelocityandPressure();
};
#endif