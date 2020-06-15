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

// initial value of resistance
double init_x = 0.0;

// population size
int const N = 10000;

int const nplasmid_max = 2;

// initial frequency of good plasmid 
// in infected individuals
double p_good_init = 0.0;
double n_plasmid_init = 1;

// density dependence
double kappa = 1;

// growth rate in absence of density dependence
double b = 1;

// number of timesteps that the simulation should run
int max_time = 10;

// individual
struct Individual
{
    // resistance
    double x;
    
    // list of good vs bad plasmids
    // note that we allow potentially for 
    // larger number of plasmids / cell than just 2
    bool plasmid_good[nplasmid_max];

    // number of plasmids carried by individual
    int nplasmids;

    // actual plasmid phenotype is then the fraction
    // of good plasmids (i.e., additive effects)
    double fraction_good;
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
    // auxiliary variable to keep track
    // of the fraction of good plasmids in each infected
    int n_good;

    // initialize susceptibles
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        Susceptible[S_idx].x = init_x;
        Susceptible[S_idx].nplasmids = 0;
    }// for (int S_idx = 0; S_idx < Ns; ++S_idx)

    // initialize infected individuals
    for (int I_idx = 0; I_idx < Ni; ++I_idx)
    {
        Infected[I_idx].x = init_x;
        Infected[I_idx].nplasmids = 1;

        // initialize fraction of good plasmids
        fraction_good = 0.0;

        for (int plasmid_idx = 0; plasmid_idx < n_plasmid_init; ++plasmid_idx)
        {
            if (uniform(rng_r) < p_good_init)
            {
                ++n_good;
                // determine whether individual carries good / bad plasmid
                Infected[I_idx].plasmid_good[plasmid_idx] = true;
            }
            else
            {
                Infected[I_idx].plasmid_good[plasmid_idx] = false;
            }

            Infected[I_idx].fraction_good = (double)n_good/n_plasmid_init;
        }
    }
}//end void 

// initialize the parameters through the command line
void init_arguments(int argc, char ** argv)
{
    max_time = atof(argv[1]);
    p_good_init = atof(argv[2]);
    kappa = atof(argv[3]);
    b = atof(argv[4]);
    base_name = argv[5];

}//end init_arguments()

void write_parameters(std::ofstream &data_file)
{
    data_file << std::endl << std::endl
        << "p_good_init" << ";" << p_good_init << std::endl
        << "N" << ";" << N << std::endl
        << "max_time" << ";" << max_time << std::endl
        << "b" << ";" << b << std::endl
        << "kappa" << ";" << kappa << std::endl
        << "n_plasmid_init" << ";" << p_good_init << std::endl;

} // end write_parameters()

// birth event of an individual
void birth(Individual &parent)
{
}

// write headers to the datafile
void write_data_headers(std::ofstream &data_file)
{
    std::vector<double> 
} // end write_data_headers()

// setup a distribution of events and choose
// which event to do
void event_chooser(int const time_step)
{
    // cumulative birth probability
    double cumul_birth = 0.0;

    // initialize a vector that contains all the weights 
    // re an individual's birth
    std::vector <double> birth_weights;

    // go through all susceptibles and calculate birth probs
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        birth_weights.push_back(prob_birth_susceptible);
    }

    // birth of a susceptible host (from another one) 
    // is easy as we do not need to
    // track any plasmids.

    // now go through the infected hosts and determine
    // their plasmid phenotype 
    for (int inf_idx = 0; inf_idx < Ni; ++inf_idx)
    {
        fecundity(Infected[inf_idx].fraction_good) * b
    }
}

int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    
    // initialize file to write data to
    std::string file_name = base_name + ".csv";
    std::ofstream data_file{file_name};
    write_data_headers(data_file);

    // initialize the population
    init_pop();

    // numerically iterate the simulation
    for (int time_idx = 0; time_idx < max_time; ++time_idx)
    {
        event_chooser(time_idx);
    }
    
    write_parameters(data_file);
}



