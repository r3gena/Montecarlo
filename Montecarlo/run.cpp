#include "Sistema_Montecarlo.h"

int Total_pasos=5000000;
int num_procesos=10;
int pasos_proceso=Total_pasos/num_procesos;

std::string energy_name;
std::string process_name;
std::string results_name="../results.txt";

int main()
{
  Sistema_MC sistema_mc("../estabilizado.txt");

  for(int i=0; i<num_procesos; i++)
  {
    energy_name="../Energy/energy_";
    process_name="Process ";

    energy_name.append(std::to_string(i+1));
    process_name.append(std::to_string(i+1));

    energy_name.append(".txt");

    sistema_mc.metropolis(pasos_proceso, energy_name, process_name, results_name);
  }
  return 0;
}
