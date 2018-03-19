#pragma once
#include <vector>

class Timing
{
public:
	int Te; //Time in hours
	int ts; //start (s)
	int te; //end (s)
	double dt; //sizes of time steps
	std::vector<double> T; //moments of time
 
	Timing();
	void SetTimeGrid(int M);
	~Timing();
};

