#pragma once
#include <vector>
#include <cmath>
#include "Fluid.h"
#include "Reservoir.h"
#include "Timing.h"
#include "BC.h"

struct Boundary
{
	int type;
	std::vector<double> bound;
};

struct ThreeDiag
{
	std::vector<double> Main;
	std::vector<double> Upper;
	std::vector<double> Lower;
};

class Solver
{
public:
	double koeff; //регулирует тип разностной схемы (0 - неявная,
	//1 - явная, 1/2 - Кранка - Николсона
	int N; //size of the grid, вспомогат.
	std::vector<double> Chi; //Piezoconductivity
	std::vector<double> KappaPlus; //Flow to the next cell
	std::vector<double> KappaMinus; //Flow to the previous cell
	std::vector<double> P; //vector of relevant pressures
	std::vector<double> Ppr; //vector of pressures on previous time step
	std::vector<double> T; //vector of relevant temperatures
	std::vector<double> Tpr; //vector of temperatures on previous time step
	Boundary Bound;
	ThreeDiag LHS;
	std::vector< std::vector<double> > B;
	std::vector<double> RHS; //Right - handed side of the equation


	Solver(Reservoir *Res, Fluid *Fl, Timing *Tm, BC *b);
	void SetPiezo(Reservoir *Res, Fluid *Fl);
	void ApprMassFlows(Reservoir *Res);
	//void FindPressure(Reservoir *Res, Fluid *Fl, Timing *T, Boundary *b);
	//std::vector<double> TDMA(std::vector<double>& D1, std::vector<double>& D2, std::vector<double>& D3, std::vector<double>& RHS);
	std::vector<double> TDMA(std::vector<double>& RHS);
	std::vector<double> Multiply(std::vector<std::vector<double>>& B, std::vector<double>& P);
	void SetLHS(Reservoir *Res, BC *b, Timing *Tm);
	void SetB(Reservoir *Res, Timing *Tm);
	~Solver();
};

