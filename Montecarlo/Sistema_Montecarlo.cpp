#include "Sistema_Montecarlo.h"

Sistema_MC::~Sistema_MC()
{
	delete[] rx;
	delete[] ry;
	delete[] rz;
}

Sistema_MC::Sistema_MC(std::string filename)
{
	//abro el archivo en cuastion
	std::ifstream inputfile(filename, std::ios::in);

	inputfile.close();
	//inicializo los generadores de numeros aleatorios
	this.int_distribution = std::uniform_int_distribution<int>(0, this.N - 1);
	this.double_distribution = std::uniform_real_distribution<double>(0, 1);
	//Para generar las configuraciones ya esta la clase sistema, esta abre lo que le mandes

	//llamo al potencial para calcular la energia correspondiente a l estado inicial
}

void Sistema_MC::metropolis(int N_steps, std::string binary_name, std::string energy_name)
{
	double old_E, new_E, deltaE, old_dphi, new_dphi, delta_dphi, old_d2phi, new_d2phi, delta_d2phi;//energias de una sola particula y variacion de energia
	int n_move;
	double xmovement, ymovement, zmovement;

	std::ofstream bin__output(binary_name, std::ios::binary);
	std::ofstream energy_output(energy_name, std::ios::out);

	for(int i=0;i<N_steps;i++)
	{
	//genero aleatoriamente que particula se mueve y cuanto
		n_move = this.int_distribution(this.generator);
		xmovement = this.double_distribution(this.generator) * this->L / 100;
		ymovement = this.double_distribution(this.generator) * this->L / 100;
		zmovement = this.double_distribution(this.generator) * this->L / 100;

		//calculamos la diferencia de energia al mover la particula con el potencial de una particula
		LJ_particula(this.N, n_move, this.rc, this.rc2, this.L, this.V, &old_E, &old_dphi, &old_d2phi, rx, ry, rz);
		//movemos la particula
		rx[n_move]+=xmovement;
		ry[n_move]+=ymovement;
		rz[n_move]+=zmovement;

		LJ_particula(this.N, n_move, this.rc, this.rc2, this.L, this.V, &new_E, &new_dphi, &new_d2phi, rx, ry, rz);
		deltaE = new_E - old_E;
		delta_dphi = new_dphi - old_dphi;
		delta_d2phi = new_d2phi - old_d2phi;

		//ahora evaluamos si el estado se descarta o se admite
		if (deltaE < 0 || this.double_distribution(this.generator) <= (exp(-deltaE / this.T)))
		{
			//Se admite, solo hay que actualizar las energias
			this.Ep +=deltaE;
			this.dphi += delta_dphi;
			this.d2phi += delta_d2phi;
		}
		else
		{
			//se rechaza por lo cual deshacemos el movimiento, las energias permanecen igual
			rx[n_move]-=xmovement;
			ry[n_move]-=ymovement;
			rz[n_move]-=zmovement;
		}
		//realizar el guardado
		//binarios
		bin__output.write((char*)rx, this.N*sizeof(double));
		bin__output.write((char*)ry, this.N*sizeof(double));
		bin__output.write((char*)rz, this.N*sizeof(double));
		//energias
		energy_output << this.Ep << " " << this.dphi << " " << this.d2phi << std::endl;
	}
	//cerramos los archivos
	bin__output.close()
	energy_output.close()
}
