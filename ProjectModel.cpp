// ProjectModel.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
#include <vector>
#include <string>
#include "Reservoir.h"
#include "Fluid.h"
#include "Solver.h"
#include <ctime>
#include "Timing.h"

using namespace std;

int main()
{
	unsigned int start_time = clock();

	//Time
	
	int const M = 200; //number of time steps
	Timing *Tm = new Timing();
	Tm->SetTimeGrid(M);
	
	Reservoir *Res = new Reservoir();
	double dr1 = 0.01;
	double dr0 = 1;

	int N1 = 50;
	int N2 = 35;
	
	dr0 = 0.1;
	dr1 = 10;
	double step_k = 1.5;
	Res->Grid_simple(dr1, dr0);
	//Res->Grid_log_reg(N1, N2);
	//Res->Grid_log_accel(0.01, 1.1, 10);
	//Res->Grid_log_accel(0.01, 1.1, 1.3, 10);
	//Res->Grid_simple_accel(dr0, step_k, dr1);
	//Res->Grid_simple_accel(0.01, 1.1, 1.4, 10);

	//system("pause");

	int N = Res->Rplus.size();

	//FLUID
	Fluid *F = new Fluid();

	//TYPE OF SCHEME
	/* IF EXPLISIT SCHEME IS USED, koeff = 1,
	IF IMPLICIT, koeff = 0;
	CRANK - NICOLSON SCHEME - koeff = 1/2.
	*/
	double koeff = 0;
	
	string fname;
	/*
	fname = "Pw(t),";
	fname += "M=";
	fname += to_string(M);
	fname += ",S=";
	fname += to_string((int)Res->s);
	fname += ",T=";
	fname += to_string(Tm->Te);
	fname += ".csv";
	*/
	fname = "AccelGrid.csv";

	ofstream fout(fname);

	fout << "Time" << ',' << "P(n)";
	fout << "\n";

	string whole_P;
	whole_P = "Whole_pressure ";
	whole_P += to_string(dr0);
	//whole_P += to_string(N1);
	whole_P += ", accel = ";
	whole_P += to_string(step_k);
	whole_P += ", ";
	whole_P += to_string(dr1);
	//whole_P += to_string(N2);
	whole_P += ".csv";

	ofstream p_out(whole_P);
	p_out << "Time (s)" << "," << Res->rw<< ",";

	for (int i = 0; i < Res->R.size(); i++) {
		p_out << Res->R[i] << ",";
	}
	p_out << Res->Rplus[Res->Rplus.size() - 1];
	p_out << "\n";

	//BOUNDARY CONDITIONS
	BC *bound = new BC(2, 1, 3 * pow(10, 7), N, Res, F);

	//PIEZOCONDUCTIVITY
	Solver *Solve = new Solver(Res, F, Tm, bound);

	/*Types of BC supported:
	bType = 1 - constant pressure on the well
	bType = 2 - constant production rate
	*/
	int bType = 2; //flag for type of BC

	double Pe = 3 * pow(10, 7);  //reservoir pressure (Pa) - used as BC in both bTypes

	//Vector of boundary conditions
	vector<double> Boundary(N + 2);
	for (int i = 1; i < (Boundary.size() - 1); i++) {
		Boundary[i] = 0;
	}
	Boundary[N + 1] = Pe;

	//Matrix of pressures: array of vectors

	vector < double > P[M + 1];

	for (int i = 0; i < M + 1; i++) {
		P[i].resize(N + 2);
	}

	for (int i = 0; i < N + 2; i++) {
		P[0][i] = Pe;  // at t=0 pressure along the reservoir is equal to Pe
		for (int j = 1; j < M + 1; j++) {
			P[j][i] = 0.; // those pressures we don't know yet
			P[j][N + 1] = Pe; // constant pressure at the edge
		}
	}

	double Pw, Q;

	// Differences raised by different boundary conditions types

	switch (bType) {
	case 1: {
		Pw = 1.2 * pow(10, 7);  //wellbore pressure (Pa) - used in bType 1 only
		Boundary[0] = Pw;  //wellbore pressure (Pa)
		for (int i = 0; i < M + 1; i++) {
			P[i][0] = Pw; // constant wellbore pressure as set
		}
		break;
	}
	case 2: {
		Q = 1; //production rate (m3/DAY!!) - used in bType 2 only
		Boundary[0] = Q * (F->mu / Res->k1 / 2 / Res->h)*log(Res -> R[0] / Res->rw) / 24 / 60 / 60 / M_PI; //production rate (m3/day)
		break;
	}
	}


	for (int i = 0; i < M + 1; i++) {
		for (int j = 0; j < N + 2; j++) {
			//cout << P[i][j] << " ";
		}
		//cout << "\n";
	}

	//RHS (Diagonal matrix B)

	vector< vector<double> > B(N + 2, vector<double>(N + 2));

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

	//for (int i = 0; i < N + 2; i++) {
	//	for (int j = 0; j < N + 2; j++) {
	//		cout << B[i][j] << " ";
	//	}
	//	cout << "\n";
	//}

	/*
	p_out << Tm->T[0] << ',';

	//p_out << RHS[0] << ",";
	for (int j = 0; j < Res->R.size(); ++j)
	{
		p_out << P[0][j] << ',';
	}
	p_out << Boundary[Boundary.size() - 1] << "\n";
	*/
	vector<double> RHS;
	//cout << "Boundary " << Boundary[0] << "\n";
	for (int t = 1; t < M+1; t++) {
		RHS = Solve->Multiply(B, P[t - 1]);
		
		

		for (int i = 0; i < RHS.size(); i++) {
			RHS[i] += Solve->Bound.bound[i];
			//cout << "RHS now" << RHS[i] << " ";
		}

		//cout << RHS[0];


		//P[t] = Solve->TDMA(LHS_Main,LHS_Upper, LHS_Lower, RHS);
		P[t] = Solve->TDMA(RHS);

		p_out << Tm->T[t] << ',';
		p_out << P[t][0] - RHS[0] << ",";
		//p_out << RHS[0] << ",";
		for (int j = 0; j < Res->R.size(); ++j)
		{
			p_out << P[t][j] << ',';
		}
		p_out <<Boundary[Boundary.size()-1]<< "\n";
		std::cout << t;

	}

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	
													  
	//THEORY

	double t_begin = pow(Res->rw, 2)*F->mu*Res -> c*Res->phi / Res->k;
	double t_end = pow(Res->re, 2)*F->mu*Res->c*Res->phi / Res->k;

	if (Res->s > 0) {
		t_end = t_end / 2;
	}
	if (Res ->s==0) {
		t_end = t_end / 4;
	}

	int num1, num2;

	for (int i = 0; i < M + 1; i++) {
		if (Tm->T[i] >= t_begin) {
			num1 = i;
			break;
		}
	}

	if (t_end > Tm->T[M]) {
		num2 = M;
	}
	else
	{
		for (int i = M; i >= num1; --i) {
			num2 = i;
			if (t_end >= Tm->T[i]) { break; };
		}
	}


	// for each row
	for (int i = num1; i < num2; ++i) {
		// for each column
		//for (int j = 0; j < data[i].size(); ++j)
		fout << Tm->T[i] << ',' << P[i][0];
		fout << "\n";
	}
	fout.close();

	//cout << fname;

	string fname1;
	fname1 = "P(log r),";
	fname1 += "M=";
	fname1 += to_string(M);
	fname1 += ",S=";
	fname1 += to_string((int)Res->s);
	fname1 += ",T=";
	fname1 += to_string(Tm->Te);
	fname1 += ".csv";

	std::cout << "\n runtime =" << search_time/1000;
	std::system("pause");
	return 0;
}

