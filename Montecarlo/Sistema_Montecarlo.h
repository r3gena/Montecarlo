#pragma once
#include "potenciales.h"
#include "progress_bar.h"//igual no se puede usar esto en este proyecto si uso cuda
#include <random>

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

	void metropolis();//es mejor declarar esta funcion entera en gpu o solo la suma de pares del potencial?
	//declarando solo el potencial en gpu es mas facil pero igual tarda mas pasar variables entre un sitio y otro que el calculo en si.

	Sistema_MC(int N, double T, double V, double rc);//no se si necesito mas variables
	~Sistema_MC();

private:
	//posiciones de las particulas
	double* rx, * ry, * rz;
	//posiciones del paso de comprobacion
	double* temp_rx, * temp_ry, * temp_rz;//igual es mejor mover y deshacer el movimiento si no se acepta que trabajar con el doble de variables
};