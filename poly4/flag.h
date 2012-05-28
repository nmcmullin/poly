#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;
#ifndef FLAG_H
#define FLAG_H
class flag
{
public:
	flag(bool inC, bool inR, bool inL, bool inB, bool inT, bool inF, bool inS);
	~flag();
	bool C,R,L,B,T,F,S;
	//Control Volume Cell? (True or False)
	//Cell to left/right/top/bottom: Control Volume cell?
	//Fluid in the Cell? (True or False)
	//Fluid Cell containing an interface or Surface (True or False)
};
typedef vector<flag> flaglist;
#endif