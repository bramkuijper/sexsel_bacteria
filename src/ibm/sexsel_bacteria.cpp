// sexual selection in bactiera

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>


// C++ random number generation unsigned int seed = get_nanoseconds();
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r{seed};
std::uniform_real_distribution<> uniform(0.0,1.0);

// file base name
std::string base_name;

// population size
int const N = 10000;

int const nplasmid_max = 2;

// number of timesteps that the simulation should run
int number_timesteps = 10;

// individual
struct Individual
{
    // resistance
    double x;
    
    // list of good vs bad plasmids
    bool plasmid_good[nplasmid_max];

    // number of plasmids carried by individual
    int nplasmids;
};

typedef Individual Population[N];

Population Susceptible;
Population Infected;


// number of infected and susceptible hosts
int Ni = 100;
int Ns = N - Ni;


// initialize population
void init_pop()
{
    // initialize susceptibles
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        Susceptible[S_idx].x = init_x;
        Susceptible[S_idx].nplasmids;
    }// for (int S_idx = 0; S_idx < Ns; ++S_idx)

    // initialize infected individuals
    for (int I_idx = 0; I_idx < Ni; ++I_idx)
    {
        Infected[I_idx].x = init_x;
        Infected[I_idx].nplasmids = 1;

        for (int plasmid_i = 0; n_plasmid_init; ++n_plasmid_init)
        {
            // randomly determine whether individual carries good / bad plasmid
            Infected[I_idx].plasmid_good[plasmid_idx] = uniform(rng_r) < 0.5;
        }

        Infected[I_idx].plasmid_good[
    }
}//end void 

int main(int argc, char **argv)
{
    init_arguments();
    init_pop();
}



