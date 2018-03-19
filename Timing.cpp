#include "stdafx.h"
#include "Timing.h"
#include <iostream>


Timing::Timing()
{
	Te = 750;
	//M = 20;
	ts = 0;
	te = Te * 60 * 60;
}

void Timing::SetTimeGrid(int M)
{
	dt = (te - ts) / M;
	T.resize(M + 2);
	for (int i = 0; i < M + 2; i++) {
		T[i] = ts + dt * i;
		//std::cout << i << "," << T[i] << "\n";
	}
	std::cout << "Time grid built\n";
}

Timing::~Timing()
{
}
