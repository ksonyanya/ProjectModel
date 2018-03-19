#include "stdafx.h"
#include "Reservoir.h"
#include <iostream>

Reservoir::Reservoir()
{	re = 300;
	rw = 0.108;
	rd = rw + 0.3;
	s = 5;
	phi = 0.2;
	c = pow(10, -9);
	h = 1;
	k = pow(10, -13) / 10;
	k1 = k / ((s / log((rd) / rw)) + 1);
	R = std::vector<double>{};
	Rplus = std::vector<double>{};
	Rminus = std::vector<double>{};
	alpha = std::vector<double>{};
}

void Reservoir::Grid_simple(double dr1, double dr0)
{
	int N1 = ceil((rd - rw) / dr1);
	dr1 = (rd - rw) / N1;
	int N2 = ceil((re - rd) / dr0);
	dr0 = (re - rd) / N2;
	int N = N1 + N2;
	std::vector<double> dr(N);

	for (int i = 0; i < N1; i++) {
		dr[i] = dr1;    //in damage zone

		for (int i = N1; i < N; i++) {
			dr[i] = dr0;    //out of damage zone
		}
	}
	FormGrid(dr);
	std::cout << dr.size() << "\n";

}

void Reservoir::Grid_log_reg(const int N1, const int N2)
{
	int N = N1 + N2;
	std::vector<double> dr(N);

	for (int i = 0; i < N1; i++) {
		dr[i] = log(rd / rw) / N1;    //in damage zone
		//std::cout << dr[i] << " ";
	}
	for (int i = N1; i < N; i++) {
			dr[i] = log(re / rd) / N2;    //out of damage zone
			//std::cout << dr[i] << " ";
	}
	FormGrid_log(dr);
	//return dr;
}

void Reservoir::Grid_log_accel(double dr0, const double k, double dr_max)
{
	std::vector<double> dr;
	dr.push_back(dr0);
	double sum_dr = dr.back();
	while (sum_dr < log(rd/rw))
	{
		dr.push_back(k*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > log(rd / rw)) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(log(rd/rw) - sum_dr);
		sum_dr += dr.back();
	}

	dr.push_back(dr.back());
	sum_dr = dr.back();
	while (sum_dr < log(re / rd))
	{
		dr.push_back(k*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > log(re / rd)) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(log(re / rd) - sum_dr);
		sum_dr += dr.back();
	}
	
	//for (int i = 0; i < dr.size(); i++) {
	//	std::cout << dr[i] << " ";
	//}
	//std::cout << "\n";
	//std::cout << sum_dr << " " << log(re/rd) << "\n";

	FormGrid_log(dr);
	std::cout << dr.size() << "\n";
}

void Reservoir::Grid_log_accel(double dr0, const double k0, const double k1, double dr_max)
{
	std::vector<double> dr;
	dr.push_back(dr0);
	double sum_dr = dr.back();
	while (sum_dr < log(rd / rw))
	{
		dr.push_back(k0*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > log(rd / rw)) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(log(rd / rw) - sum_dr);
		sum_dr += dr.back();
	}

	dr.push_back(dr.back());
	sum_dr = dr.back();
	while (sum_dr < log(re / rd))
	{
		dr.push_back(k1*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > log(re / rd)) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(log(re / rd) - sum_dr);
		sum_dr += dr.back();
	}

	//for (int i = 0; i < dr.size(); i++) {
	//	std::cout << dr[i] << " ";
	//}
	//std::cout << "\n";
	std::cout << sum_dr << " " << log(re/rd) << "\n";

	FormGrid_log(dr);
	std::cout << dr.size() << "\n";
}

void Reservoir::Grid_simple_accel(double dr0, const double k, double dr_max)
{
	std::vector<double> dr;
	dr.push_back(dr0);
	double sum_dr = dr.back();
	while (sum_dr < rd-rw)
	{
		dr.push_back(k*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > rd-rw) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(rd-rw - sum_dr);
		sum_dr += dr.back();
	}

	dr.push_back(dr.back());
	sum_dr = dr.back();
	while (sum_dr < re-rd)
	{
		dr.push_back(k*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > re-rd) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(re-rd - sum_dr);
		sum_dr += dr.back();
	}

	//for (int i = 0; i < dr.size(); i++) {
	//	std::cout << dr[i] << " ";
	//}
	//std::cout << "\n";
	//std::cout << sum_dr << " " << log(re/rd) << "\n";

	FormGrid(dr);
	std::cout << dr.size() << "\n";

}

void Reservoir::Grid_simple_accel(double dr0, const double k0, const double k1, double dr_max)
{
	std::vector<double> dr;
	dr.push_back(dr0);
	double sum_dr = dr.back();
	while (sum_dr < rd-rw)
	{
		dr.push_back(k0*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > rd-rw) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(rd-rw - sum_dr);
		sum_dr += dr.back();
	}

	dr.push_back(dr.back());
	sum_dr = dr.back();
	while (sum_dr < re-rd)
	{
		dr.push_back(k1*dr.back());
		sum_dr += dr.back();
	}

	if (sum_dr > re-rd) {
		sum_dr -= dr.back();
		dr.pop_back();
		dr.push_back(re-rd - sum_dr);
		sum_dr += dr.back();
	}

	//for (int i = 0; i < dr.size(); i++) {
	//	std::cout << dr[i] << " ";
	//}
	//std::cout << "\n";
	//std::cout << sum_dr << " " << log(re/rd) << "\n";

	FormGrid(dr);
	std::cout << dr.size() << "\n";
}

void Reservoir::FormGrid(std::vector<double> dr)
{
	int N = dr.size();
	R.resize(N);
	Rplus.resize(N);
	Rminus.resize(N);

	Rplus[0] = rw + dr[0];
	Rminus[0] = rw;

	for (int i = 1; i < Rplus.size(); i++) {
		Rplus[i] = Rplus[i - 1] + dr[i];
		Rminus[i] = Rplus[i - 1];
	}

	for (int i = 0; i < Rplus.size(); i++) {
		R[i] = Rminus[i] + (Rplus[i] - Rminus[i]) / 2;
	}

	alpha.resize(N);
	for (int i = 0; i < N; i++) {
		alpha[i] = ((pow(Rplus[i], 2.0) - pow(Rminus[i], 2.0)) / 2);
		//std::cout << R[i] << " ";
	}

	//std::cout << "\n";
	//std::cout << dr.size() << "\n";
}

void Reservoir::FormGrid_log(std::vector<double> dr)
{
	int N = dr.size();
	R.resize(N);
	Rplus.resize(N);
	Rminus.resize(N);

	Rplus[0] = log(rw) + dr[0];
	Rminus[0] = log(rw);

	//std::cout << Rminus[0] << " " << Rplus[0] << "\n";
	//system("pause");

	for (int i = 1; i < Rplus.size(); i++) {
		Rplus[i] = Rplus[i - 1] + dr[i];
		Rminus[i] = Rplus[i - 1];
	}

	for (int i = 0; i < Rplus.size(); i++) {
		Rplus[i] = exp(Rplus[i]);
		Rminus[i] = exp(Rminus[i]);
		//std::cout << Rminus[i] << " " << Rplus[i] << "\n";
	}

	//system("pause");
	for (int i = 0; i < Rplus.size(); i++) {
		R[i] = Rminus[i] + (Rplus[i] - Rminus[i]) / 2;
	}

	alpha.resize(N);
	for (int i = 0; i < N; i++) {
		alpha[i] = ((pow(Rplus[i], 2.0) - pow(Rminus[i], 2.0)) / 2);
		//std::cout << R[i] << " ";
	}

	//std::cout << "\n";
	//std::cout << dr.size() << "\n";
}

Reservoir::~Reservoir()
{
}
