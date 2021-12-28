#include "Sistema_Montecarlo.h"

Sistema_MC::sistema("estabilizado.txt");

std::string energy_name="Energy/energy_test.txt";
std::string binary_name="Binary/binary_test.txt";
int N_pasos=100000

sistema.metropolis(N_pasos, binary_name, energy_name, "prueba");
