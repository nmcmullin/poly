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
#ifndef FreeSurfaceCalc_H
#define FreeSurfaceCalc_H
class FreeSurfaceCalc : virtual public CFD
{
public:
	FreeSurfaceCalc(int inncols, int innrows, float indx, float indy, float inRe, float inrho, float indynvis, bool inSMACbool,int innpart,bool inSLICbool,bool inFLAPbool,bool inYVOFbool, float invoftol);
	virtual ~FreeSurfaceCalc();
protected:
	virtual void RunFreeSurfaceCalc();
	virtual void UpdateFluidFlag();
	virtual void PrintFreeSurfaceData(int &myMval);

private:
	void UpdateParticleLocations();
	bool SetFluidFlagTrue(int &passI, int &passJ);
	bool SetFluidFlagFalse(int &passI, int &passJ);
	void SetEmptyCellstoZero();
	void YVOFInterfaceConstantCalc(int &I, int &J, float &ALPHA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcXsweepA(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcXsweepBIIrevIII(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcXsweepC(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcXsweepDIIrevIII(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcYsweepA(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcYsweepBIIrevIII(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcYsweepC(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void FluxCalcYsweepDIIrevIII(int &I, int &J,float &ALPHA, float &BETA, float &ST, float &SR, float &SB, float &SL);
	void RunSLIC();
	void SLICXsweep(int &I, int &J);
	void SLICYsweep(int &I, int &J);
	void RunFLAP();
	void FLAPXsweep(int &I, int &J);
	void FLAPYsweep(int &I, int &J);
	void RunYVOF();
	void RunSMAC();
	void PeriodicVolumeFractionBC();
	void VolumeFractionBC();
	void ImposePeriodicFluxBC();
	int const npart;
	float voftol;//Tolerance for deciding that C is close enough to 1 or zero, for it to be considered ==1 or ==0
	bool initialdirectionbool;//boolean variable keeping track of sweep direction balance
	vector<float> Fr,Fl,Ft,Fb;//Fluxes out of each cell's right edge, left edge, top edge and bottom edge
};
#endif
