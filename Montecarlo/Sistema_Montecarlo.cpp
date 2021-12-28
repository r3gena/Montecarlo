#include "Sistema_Montecarlo.h"

Sistema_MC::~Sistema_MC()
{
	delete[] rx;
	delete[] ry;
	delete[] rz;
}

Sistema_MC::Sistema_MC(int N, double T, double V, double rc)
{
	//inicializo los generadores de numeros aleatorios
	// std::default_random_engine generator;
	this->int_distribution = std::uniform_int_distribution<int>(0, this->N - 1);
	this->double_distribution = std::uniform_real_distribution<double>(0, 1);
	//primero hay que colocar las particulas en una configuracion inicial

	//llamo al potencial para calcular la energia correspondiente a l estado inicial
}

void Sistema_MC::metropolis()
{
	double old_E, new_E, deltaE, old_dphi, new_dphi, delta_dphi, old_d2phi, new_d2phi, delta_d2phi;//energias de una sola particula y variacion de energia

	//genero aleatoriamente que particula se mueve y cuanto
	int n_move = this->int_distribution(this->generator);
	double xmovement = this->double_distribution(this->generator) * this->L / 100;
	double ymovement = this->double_distribution(this->generator) * this->L / 100;
	double zmovement = this->double_distribution(this->generator) * this->L / 100;

	//calculamos la diferencia de energia al mover la particula con el potencial de una particula
	LJ_particula(N, n_move, rc, rc2, L, V, &old_E, &old_dphi, &old_d2phi, rx, ry, rz);
	//movemos la particula
	rx[n_move]+=xmovement;
	ry[n_move]+=ymovement;
	rz[n_move]+=zmovement;

	LJ_particula(N, n_move, rc, rc2, L, V, &new_E, &new_dphi, &new_d2phi, rx, ry, rz);
	deltaE = new_E - old_E;
	delta_dphi = new_dphi - old_dphi;
	delta_d2phi = new_d2phi - old_d2phi;

	//ahora evaluamos si el estado se descarta o se admite
	if (deltaE < 0 || this->double_distribution(this->generator) <= (exp(-deltaE / this->T)))
	{
		//Se admite

		this->Ep +=deltaE;
		this->dphi += delta_dphi;
		this->d2phi += delta_d2phi;
		//realizar el guardado
	}
	else
	{
		//se rechaza por lo cual deshacemos el movimiento, las energias permanecen igual
		rx[n_move]-=xmovement;
		ry[n_move]-=ymovement;
		rz[n_move]-=zmovement;
		//guardar estado actual de nuevo
	}
}
