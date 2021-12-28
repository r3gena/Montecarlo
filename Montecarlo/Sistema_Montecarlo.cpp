#include "Sistema_Montecarlo.h"

Sistema_MC::~Sistema_MC()
{
	delete[] rx;
	delete[] ry;
	delete[] rz;
}

Sistema_MC::Sistema_MC(std::string filename)
{
	//abro el archivo
	std::string bin_filename;
	std::ifstream inputfile(filename, std::ios::in);
	inputfile >> this->N >> this->rho >> this->L >> this->rc;
  inputfile >> this->T;
  inputfile >> bin_filename;
	inputfile.close();

	this->V=L*L*L;
	this->rc2=rc*rc;
	//inicializo los generadores de numeros aleatorios
	this->int_distribution = std::uniform_int_distribution<int>(0, N - 1);
	this->double_distribution = std::uniform_real_distribution<double>(0, 1);

	//inicializo las posiciones
	this->rx = new double[N];
	this->ry = new double[N];
	this->rz = new double[N];

	std::ifstream binfile(bin_filename, std::ios::binary);
	binfile.read((char*)rx, N*sizeof(double));
	binfile.read((char*)ry, N*sizeof(double));
	binfile.read((char*)rz, N*sizeof(double));
	binfile.close();
	//llamo al potencial para calcular la energia correspondiente al estado inicial
	lennard_jones(N, rc, rc2, L, V, &Ep, &dphi, &d2phi, rx, ry, rz);
}

void Sistema_MC::metropolis(int N_steps, std::string energy_name, std::string barname, std::string results_name)
{
	double old_E, new_E, deltaE, old_dphi, new_dphi, delta_dphi, old_d2phi, new_d2phi, delta_d2phi;//energias de una sola particula y variacion de energia
	int n_move;
	double xmovement, ymovement, zmovement;

	std::ofstream energy_output(energy_name, std::ios::out);

	//creo la barra de progreso
	progress_bar bar(barname, N_steps);

	//medias que queremos calcular:
	double mphi=0;//potencial
	double mdphi=0, md2phi=0;//primera y segunda derivada
	double mphi2=0, mdphi2;//potencial y primera derivada al cuadrado
	double mphidphi=0;//potencial por su primera derivada

	for(int i=0;i<N_steps;i++)
	{
	//genero aleatoriamente que particula se mueve y cuanto
		n_move = int_distribution(generator);
		xmovement = double_distribution(generator) * L / 100;
		ymovement = double_distribution(generator) * L / 100;
		zmovement = double_distribution(generator) * L / 100;

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
		bool test1=deltaE<0;
		bool test2=double_distribution(generator)<=(exp(-deltaE/T));
		bool test_global= test1||test2;
		if (test_global)
		{
			//Se admite, solo hay que actualizar las energias
			Ep +=deltaE;
			dphi += delta_dphi;
			d2phi += delta_d2phi;
		}
		else
		{
			//se rechaza por lo cual deshacemos el movimiento, las energias permanecen igual
			rx[n_move]-=xmovement;
			ry[n_move]-=ymovement;
			rz[n_move]-=zmovement;
		}
		//realizar el guardado
		energy_output << Ep << " " << dphi << " " << d2phi << std::endl;

		//añadimos los valores aceptados a las medias
		mphi=mphi+Ep;
		mdphi=mdphi+dphi;
		mphi2=mphi2+Ep*Ep;
		mdphi2=mdphi2+dphi*dphi;
		md2phi=md2phi+d2phi;
		mphidphi=mphidphi+Ep*dphi;

		bar.update(i);
	}
	//cerramos el archivo
	energy_output.close();

	//dividimos las medias por el numero de pasos
	mphi=mphi/N_steps;
	mdphi=mdphi/N_steps;
	mphi2=mphi2/N_steps;
	mdphi2=mdphi2/N_steps;
	md2phi=md2phi/N_steps;
	mphidphi=mphidphi/N_steps;

	//calculamos los coeficientes y los guardamos al archivo correspondiente
	double p=(N*T/V)-mdphi;
	double Cv=((3*N-3)/2)+(mphi2-mphi*mphi)/(T*T);
	double gamma=(V/Cv)*((N/V)+(mphi*mdphi-mphidphi)/(T*T));
	double inv_kT=(N*T/V)+V*md2phi-(V/T)*(mdphi2-mdphi*mdphi);

	std::ofstream results_file(results_name, std::ios::app);
	results_file << p << " " << Cv << " " << gamma << " " << 1/inv_kT << std::endl;
	results_file.close();
}
