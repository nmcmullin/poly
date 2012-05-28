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
#include "CFD.h"
#define PI 3.14159
using namespace std;
#ifndef NavierStokesCalc_H
#define NavierStokesCalc_H
class NavierStokesCalc : virtual public CFD
{
public:
	NavierStokesCalc(int inncols, int innrows, float indx, float indy, float inRe, float inrho, float indynvis, bool inSMACbool, bool inSLICbool, bool inFLAPbool, bool inYVOFbool, float intmax, int inmaxNumIt, float inveleps, float inpeps, float ingx, float ingy, float inw, float intsf, float ingsf, bool inGSbool, bool inCGbool);
	virtual ~NavierStokesCalc();
protected:
	virtual void RunNavierStokes();
	virtual	void SetPreliminaryGuessValues();
	virtual bool UpdateVelocity();
	float const tmax;
private:
	void DynamicAssignDTGAM();
	void CalculatePreliminaryGuess();
	void Fcalc(int &i, int &j);
	void Gcalc(int &i, int &j);
	void ImposevelocityBC();
	void ImposeSpecialVelocityBC();
	void ImposePeriodicVelocityBC();
	void ImposePeriodicPressureBC();
	void ImposepressureBC();
	void ImposePeriodicpressureBC();
	void CalcRHS();
	bool PressureResidualCheck();
	void SORpressure();
	void InitialConditionsForCG();
	void InitialConditionsForCG2();
	void CGpressure();
	float const gx,gy,w,veleps,peps,tsf,gsf;
	int const maxNumIt;
	float gam,rsum,newrsum,sumpd,alpha,beta;
	vector<float> F,G,rhs,pd,r;
	bool GSbool,CGbool;
};
#endif