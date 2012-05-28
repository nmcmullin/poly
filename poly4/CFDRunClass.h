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
#include "NavierStokesCalc.h"
#include "FreeSurfaceCalc.h"
#define PI 3.14159
using namespace std;
#ifndef CFDRunClass_H
#define CFDRunClass_H
class CFDRunClass : public NavierStokesCalc, public FreeSurfaceCalc
{
public:
	CFDRunClass(int inncols,int innrows, float indx, float indy, float inRe, float inrho, float indynvis, float intmax, float inmaxNumIt, float inveleps, float inpeps, float ingx, float ingy, float inw, float intsf, float ingsf, bool inGSbool, bool inCGbool, bool inSMACbool, int innpart, bool inSLICbool, bool inFLAPbool, bool inYVOFbool, float invoftol);
	virtual ~CFDRunClass();
	void Run();
};

#endif