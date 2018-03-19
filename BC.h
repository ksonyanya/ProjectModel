#pragma once
#include <vector>
#include "Reservoir.h"
#include "Fluid.h"

class BC
{
public:
	int type;
	std::vector<double> bound;

	BC(int type, double Cond_Right, double Pe, int N);
	BC(int type, double Q, double Pe, int N, Reservoir *Res, Fluid *Fl);
	void PrintBC();
	~BC();
};

