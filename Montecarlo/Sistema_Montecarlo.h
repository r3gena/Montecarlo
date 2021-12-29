#pragma once
#include "potenciales.h"
#include "progress_bar.h"//igual no se puede usar esto en este proyecto si uso cuda
#include <random>
#include <ostream>//por algun motivo tienen nombres distintos aqui????
#include <istream>
#include <string>

class Sistema_MC
{
public:
	int N;
	double Ep, dphi, d2phi, T;
	double rc, rc2, L, V, rho, dr;

	//elementos para los numeros aleatorios
	std::default_random_engine generator;
	std::uniform_int_distribution<int> int_distribution;
	std::uniform_real_distribution<double> double_distribution;

	void metropolis(int N_steps, std::string energy_name, std::string barname, std::string results_name);

	Sistema_MC(std::string filename);
	~Sistema_MC();

private:
	//posiciones de las particulas
	double* rx, * ry, * rz;
};
