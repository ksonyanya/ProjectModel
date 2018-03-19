#include "stdafx.h"
#define _USE_MATH_DEFINES
#include "BC.h"
#include <iostream>

BC::BC(int type, double Cond_Right, double Pe, int N)
{
	//N = KappaMinus.size();
	bound.resize(N + 2);

	//setting all elements of boundary conditions as 0 except for the last one
	//(which is the countour pressure) and the first one (which is the condition
	//on the wellbore
	for (int i = 1; i < (bound.size() - 1); i++) {
		bound[i] = 0;
	}
	bound[N + 1] = Pe;
	bound[0] = Cond_Right;

	this->type = type;
	//std::cout << "\n BC type is " << Bound.type << "\n";
	//for (int i = 0; i < (Bound.bound.size()); i++) {
	//	std::cout << Bound.bound[i] << " ";
	//}
	//std::cout << "\n";
}
BC::BC(int type, double Q, double Pe, int N, Reservoir *Res, Fluid *Fl) {
	//N = KappaMinus.size();
	bound.resize(N + 2);

	//setting all elements of boundary conditions as 0 except for the last one
	//(which is the countour pressure) and the first one (which is the condition
	//on the wellbore
	for (int i = 1; i < (bound.size() - 1); i++) {
		bound[i] = 0;
	}
	bound[N + 1] = Pe;
	bound[0] = Q * (Fl->mu / Res->k1 / 2 / Res->h)*log(Res->R[0] / Res->rw) / 24 / 60 / 60 / M_PI; //production rate (m3/day)

	this->type = type;
	//PrintBC();
}

void BC::PrintBC() {
	std::cout << "\n BC type is " << type << "\n";
	for (int i = 0; i < (bound.size()); i++) {
		std::cout << bound[i] << " ";
	}
	std::cout << "\n";
}

BC::~BC()
{
}
