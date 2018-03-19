#pragma once
#include <vector>
#include <cmath>

class Reservoir
{
public:
	double re; 
	double rw; //well radius (m)
	double rd; //damage-zone radius (m)
	double s; //skin factor (1)
	double k; //permeability within the reservoir (m^2)
	double k1; //permeability within the damage zone (m^2)
	double phi; //porosity (1)
	double c; //compressibility of the reservoir (1/Pa)
	double h; //thickness (m)
	std::vector<double> R; //array of nodes
	std::vector<double> Rplus, Rminus; //arrays of left and right edges of cells
	std::vector<double> alpha; //vector used in FVM
	//std::vector<double> Chi;

	Reservoir();
	void Grid_simple(double dr1, double dr0);
	void Grid_log_reg(const int N1, const int N2);
	void Grid_log_accel(double dr0, const double k, double dr_max);
	void Grid_log_accel(double dr0, const double k0, const double k1, double dr_max);
	void Grid_simple_accel(double dr0, const double k, double dr_max);
	//dr0 (м) - начальный шаг
	//k0 - коэффициент разгона внутри damage-зоны: dr[i+1] = k*dr[i]
	//k1 - коэффициент разгона вне damage-зоны.
	//dr_max - максимальный возможный шаг: после того, как dr[N] = dr_max, все dr[k], k>N dr[k] = dr_max
	void Grid_simple_accel(double dr0, const double k0, const double k1, double dr_max);
	void FormGrid(std::vector<double> dr);
	void FormGrid_log(std::vector<double> dr);
	//std::vector<double> TDMA(std::vector<double>& D1, std::vector<double>& D2, std::vector<double>& D3, std::vector<double>& RHS);
	//std::vector<double> Multiply(std::vector<std::vector<double>>& B, std::vector<double>& P);
	


	//std::vector <double> GetVector(int n);
	~Reservoir();
};

