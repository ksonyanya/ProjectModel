#include "stdafx.h"
#define _USE_MATH_DEFINES
#include "Solver.h"
#include <iostream>
#include <cmath>



Solver::Solver(Reservoir *Res, Fluid *Fl, Timing *Tm, BC *b)
{
	SetPiezo(Res, Fl);
	N = Res->Rplus.size();
	std::cout << N;
	ApprMassFlows(Res);
	Bound.bound = b->bound;
	Bound.type = b->type;
	if (b->type == 1) 
	{
		Ppr[0] = b->bound[0];
		P[0] = b->bound[0];
	}
	SetLHS(Res, b, Tm);
	SetB(Res, Tm);
	RHS = Multiply(B, Ppr);
	for (int i = 0; i < RHS.size(); i++) {
		RHS[i] += Bound.bound[i];
		std::cout << RHS[i] << " ";
	}
}

void Solver::SetPiezo(Reservoir *Res, Fluid *Fl)
{
	Chi.resize(Res->Rplus.size());
	for (int i = 0; i < Chi.size(); i++) {
		if (Res->R[i] < Res->rd) {
			Chi[i] = Res->k1 / (Fl->mu*Res->c*Res->phi);
		}
		else
		{
			Chi[i] = Res->k / (Fl->mu*Res->c*Res->phi);
		}
		//std::cout << Res -> Rplus[i] << " " << Chi[i] << "\n";
	}
	std::cout << "Piezoncductivity set\n";
	
}

void Solver::ApprMassFlows(Reservoir *Res) 
{
	int N = Res->Rplus.size();
	KappaPlus.resize(N);
	
	for (int j = 0; j < N - 1; j++) {
		KappaPlus[j] = 1 / ((1 / Chi[j])*log(Res->Rplus[j] / Res->R[j]) + (1 / Chi[j + 1])*log(Res->R[j + 1] / Res->Rplus[j]));
		//std::cout << KappaPlus[j] << " ";
	}
	KappaPlus[N - 1] = 1 / ((1 / Chi[N - 1])*log(Res->Rplus[N - 1] / Res->R[N - 1]));

	KappaMinus.resize(Res->Rplus.size());
	for (int j = 1; j < N; j++) {
		KappaMinus[j] = 1 / ((1 / Chi[j - 1])*log(Res->Rminus[j] / Res->R[j - 1]) + (1 / Chi[j])*log(Res->R[j] / Res->Rminus[j]));

	}
	KappaMinus[0] = 1 / ((1 / Chi[0])*log(Res->R[0] / Res->Rminus[0]));

	std::cout << "Flows approximated\n";
}

void Solver::SetLHS(Reservoir *Res, BC *b, Timing *Tm) {
	//std::cout << N << " ";
	LHS.Main.resize(N+2);
	LHS.Upper.resize(N+1);
	LHS.Lower.resize(N+1);

	for (int i = 0; i < N; i++) {
		LHS.Main[i + 1] = (Res->alpha[i] / Tm->dt) + (1 - koeff)*KappaPlus[i] + KappaMinus[i] * (1 - koeff);
	}

	for (int i = 0; i < N; i++) {
		LHS.Upper[i + 1] = -KappaPlus[i] * (1 - koeff);
	}

	for (int i = 0; i < N; i++) {
		LHS.Lower[i] = -KappaMinus[i] * (1 - koeff);
	}

	LHS.Main[N + 1] = 1.;
	LHS.Lower[N] = 0.;

	switch (b->type)
	{
	case 1: {
		LHS.Main[0] = 1.;
		LHS.Upper[0] = 0.;
		break;
	}
	case 2: {
		LHS.Main[0] = -1.;
		LHS.Upper[0] = 1.;
		break;
	}

	}
	std::cout << "LHS is set\n";
	//for (int i = 0; i < LHS.Main.size(); i++) {
	//	std::cout << LHS.Main[i] << " ";
	//}
}

void Solver::SetB(Reservoir * Res, Timing * Tm)
{
	B.resize(N + 2);
	for (int i = 0; i < B.size(); i++) {
		B[i].resize(N + 2);
	}

	for (int i = 1; i < N + 1; i++) {
		for (int j = 0; j < N + 2; j++) {
			if (i == j) {
				B[i][j] = Res -> alpha[i - 1] / Tm->dt;
			}
			else
			{
				B[i][j] = 0;
			}

		}
	}
}

// parameters are:
// D1 - Main diagonal of matrix A (size is equal to RHS size)
// D2 - Upper diagonal of matrix A (size is RHS size - 1)
// D3 - Lower diagonal of matrix A (size is RHS size - 1)
// Function finds the solution of the linear system of equations Ax = b,
// where matrix A is tri-diagonal.
// Algorythm: -- forward sweep: initializing c and d coefficients via D1, D2, D3 and RHS
// -- backward sweep: obtaining the solution of the system by back substitution.
//std::vector<double> Solver::TDMA(std::vector<double>& D1, std::vector<double>& D2, std::vector<double>& D3, std::vector<double>& RHS) {
std::vector<double> Solver::TDMA(std::vector<double>& RHS) {
	std::vector<double> D1 = LHS.Main;
	std::vector<double> D2 = LHS.Upper;
	std::vector<double> D3 = LHS.Lower;

	std::vector<double> result = RHS;

	//dimensionality of the system
	int N = RHS.size();

	//coefficients
	std::vector<double> c(N);
	std::vector<double> d(N);

	//forward sweep:
	c[0] = D2[0] / D1[0];
	d[0] = RHS[0] / D1[0];

	c[N - 1] = 0;

	//cout << c[0] << "\n";
	//cout << d[0] << "\n";

	for (int i = 1; i < N - 1; i++) {
		c[i] = D2[i] / (D1[i] - D3[i - 1] * c[i - 1]);
		//cout << c[i] << "\n";
	}


	for (int i = 1; i < N; i++) {
		d[i] = (RHS[i] - D3[i - 1] * d[i - 1]) / (D1[i] - D3[i - 1] * c[i - 1]);
		//cout << d[i] << "\n";
	}

	//backward sweep
	result[N - 1] = d[N - 1];

	for (int i = N - 2; i >= 0; i--) {
		result[i] = d[i] - c[i] * result[i + 1];
	}

	return result;
}

std::vector<double> Solver::Multiply(std::vector<std::vector<double>>& B, std::vector<double>& P) {

	int M = B.size();
	int N = P.size();

	std::vector<double> result(N);

	for (int t = 1; t < M; t++) {
		for (int j = 0; j < N; j++) {
			result[j] = 0;
			for (int k = 0; k < N; k++) {
				result[j] += P[k] * B[k][j];
			}
		}
	}

	return (result);
}





Solver::~Solver()
{
}
