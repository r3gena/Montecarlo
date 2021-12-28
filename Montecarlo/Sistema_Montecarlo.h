#pragma once
#include "potenciales.h"
#include "progress_bar.h"//igual no se puede usar esto en este proyecto si uso cuda
#include <random>
#include <ofstream>
#include <ifstream>

class Sistema_MC
{
public:
	int N;
	double Ep, dphi, d2phi, T;
	double rc, rc2, L, V;

	//elementos para los numeros aleatorios
	std::default_random_engine generator;
	std::uniform_int_distribution<int> int_distribution;
	std::uniform_real_distribution<double> double_distribution;

	void metropolis();

	Sistema_MC(int N, double T, double V, double rc);//no se si necesito mas variables
	~Sistema_MC();

private:
	//posiciones de las particulas
	double* rx, * ry, * rz;
};
