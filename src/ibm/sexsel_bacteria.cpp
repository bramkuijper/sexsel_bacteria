// sexual selection in bactiera

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
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
double bmax = 1;

// resistance cost factor 
double c = 1;

// loss rate of good plasmid
double gamma_G = 0.0;
// loss rate of bad plasmid
double gamma_B = 0.0;

// infection rate of good plasmid
double psi_G = 0.0;
// infection rate of bad plasmid
double psi_B = 0.0;

// death rate susceptible
double d = 0.0;

// death rate host infected only with good plasmids
double dG = 0.0;

// death rate host infected only with bad plasmids
double dB = 0.0;

// probability of co-infection
double sigma = 0.0;

// resistance mutation rate
double mu_x = 0.0;

double sdmu_x = 0.0;

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

    // actual number of plasmids carried by individual
    int nplasmids;

    // count of the number of good plasmids
    // needs to be recalculated once we make changes to
    // plasmid composition
    int nplasmids_good;

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
        Infected[I_idx].nplasmids = n_plasmid_init;

        // initialize fraction of good plasmids
        n_good = 0;

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

        }

        Infected[I_idx].nplasmids_good = n_good;
        assert(Infected[I_idx].nplasmids_good <= Infected[I_idx].nplasmids);
        Infected[I_idx].fraction_good = (double)n_good/n_plasmid_init;
    }
}//end void 

// initialize the parameters through the command line
void init_arguments(int argc, char ** argv)
{
    // call 
    max_time = atof(argv[1]);
    p_good_init = atof(argv[2]);
    kappa = atof(argv[3]);
    bmax = atof(argv[4]);
    c = atof(argv[4]);
    gamma_G = atof(argv[5]);
    gamma_B = atof(argv[6]);
    psi_G = atof(argv[7]);
    psi_B = atof(argv[8]);
    d = atof(argv[9]);
    dG = atof(argv[10]);
    dB = atof(argv[11]);
    sigma = atof(argv[12]);
    base_name = argv[13];

}//end init_arguments()

void write_parameters(std::ofstream &data_file)
{
    data_file << std::endl << std::endl
        << "p_good_init" << ";" << p_good_init << std::endl
        << "N" << ";" << N << std::endl
        << "max_time" << ";" << max_time << std::endl
        << "bmax" << ";" << bmax << std::endl
        << "c" << ";" << c << std::endl
        << "gamma_G" << ";" << gamma_G << std::endl
        << "gamma_B" << ";" << gamma_B << std::endl
        << "psi_G" << ";" << psi_G << std::endl
        << "psi_B" << ";" << psi_B << std::endl
        << "d" << ";" << d << std::endl
        << "dG" << ";" << dG << std::endl
        << "dB" << ";" << dB << std::endl
        << "sigma" << ";" << sigma << std::endl
        << "kappa" << ";" << kappa << std::endl
        << "seed" << ";" << seed << std::endl
        << "n_plasmid_init" << ";" << p_good_init << std::endl;

} // end write_parameters()

double mutation(double const x)
{
    if (uniform(rng_r) < mu_x)
    {
        std::normal_distribution <double> distribution(0.0,sdmu_x);
        return(x + distribution(rng_r));
    }

    return(x);
}

// infection event of a susceptible with a good
// or bad plasmid
void infection_susceptible(
        int const S_idx
        ,bool plasmid_good)
{
    assert(S_idx >= 0);
    assert(S_idx < Ns);

    // check indeed that the susceptible individual
    // is not infected
    assert(Susceptible[S_idx].nplasmids == 0);

    // assign newly infected to stack of infected
    // individuals
    Infected[Ni] = Susceptible[S_idx];

    // update plasmids
    Infected[Ni].nplasmids = 1;

    // if plasmid good, var set to 1 otherwise 0
    Infected[Ni].nplasmids_good = plasmid_good;
    Infected[Ni].plasmid_good[0] = plasmid_good;
    Infected[Ni].fraction_good = plasmid_good;

    ++Ni;

    // remove old susceptible (now infected)
    // by overwriting it with the susceptible
    // at the end of the stack of susceptibles
    Susceptible[S_idx] = Susceptible[Ns - 1];
    --Ns;

    // check whether population is still within bounds
    // with density dependence this should always be the case
    assert(Ni + Ns <= N);
}


// birth event of an individual
// is going to be different with 
void birth(Individual &parent, bool parent_susceptible)
{
    Individual kid;

    // mutate the resistance allele
    kid.x = mutation(parent.x);

    std::clamp(kid.x, 0.0, 1.0);

    kid.nplasmids = 0;
    kid.nplasmids_good = 0;
    kid.fraction_good = 0;

    if (parent_susceptible)
    {
        Susceptible[Ns++] = kid;
        assert(Ns + Ni <= N);
        return;
    }

    // for now assuming strict inheritance
    // later on we might want to introduce segregational effects
    kid.nplasmids = parent.nplasmids;

    // birth of infected individual
    for (int plasmid_idx = 0; plasmid_idx < parent.nplasmids; ++plasmid_idx)
    {
        // replicate plasmid to offspring
        kid.plasmid_good[plasmid_idx] = parent.plasmid_good[plasmid_idx];
        kid.nplasmids_good += kid.plasmid_good[plasmid_idx];
    }

    kid.fraction_good = (double)kid.nplasmids_good/kid.nplasmids;

    Infected[Ni++] = kid;

    return;
}

// death of randomly chosen a susceptible individual
void death_susceptible()
{
    std::uniform_int_distribution<int> susceptible_sampler(0, Ns - 1);

    int random_susceptible = susceptible_sampler(rng_r);
    Susceptible[random_susceptible] = Susceptible[Ns - 1];
    --Ns;

    assert(Ns >= 0);
    assert(Ns + Ni <= N);
}// end death_susceptible()

// death of an infected individual at location I_idx
void death_infected(int const I_idx)
{
    Infected[I_idx] = Infected[Ni - I_idx];
    --Ni;
}

void loss_plasmid(int const I_idx, bool const plasmid_good)
{
    assert(Infected[I_idx].nplasmids > 0);

    if (plasmid_good)
    {
        assert(Infected[I_idx].nplasmids_good > 0);

        for (int plasmid_idx = 0; plasmid_idx < Infected[I_idx].nplasmids; ++plasmid_idx)
        {
            if (Infected[I_idx].plasmid_good[plasmid_idx])
            {
                Infected[I_idx].plasmid_good[plasmid_idx] = 
                    Infected[I_idx].plasmid_good[Infected[I_idx].nplasmids - 1];

                --Infected[I_idx].nplasmids;
                --Infected[I_idx].nplasmids_good;
                break;
            }
        }

        if (Infected[I_idx].nplasmids == 0)
        {
            Susceptible[Ns++] = Infected[I_idx];

            //
            Infected[I_idx] = Infected[Ni - 1];
            --Ni;
        }
        else
        {
            Infected[I_idx].fraction_good  = (double) Infected[I_idx].nplasmids_good / 
                    Infected[I_idx].nplasmids;
        }

        return;
    } // end if (plasmid_good)


    // from here only bad plasmid loss
    //
    //
    assert(Infected[I_idx].nplasmids - Infected[I_idx].nplasmids_good > 0);

    for (int plasmid_idx = 0; plasmid_idx < Infected[I_idx].nplasmids; ++plasmid_idx)
    {
        if (!Infected[I_idx].plasmid_good[plasmid_idx])
        {
            Infected[I_idx].plasmid_good[plasmid_idx] = Infected[I_idx].plasmid_good[Infected[I_idx].nplasmids - 1];
            --Infected[I_idx].nplasmids;
            break;
        }
    } // end for

    if (Infected[I_idx].nplasmids == 0)
    {
        Susceptible[Ns++] = Infected[I_idx];

        // delete infected individual
        Infected[I_idx] = Infected[Ni - 1];
        --Ni;
    }
    else
    {
        Infected[I_idx].fraction_good  = (double) Infected[I_idx].nplasmids_good / 
                Infected[I_idx].nplasmids;
    }
    return;
}// end loss_plasmid()

// write headers to the datafile
void write_data_headers(std::ofstream &data_file)
{
    data_file << "time;mean_resistance;Ns;Ni;mean_nplasmid;" << std::endl;
} // end write_data_headers()

// fecundity function accounting for plasmid behaviour
// (good vs bad)
double F(double const fraction_good)
{
    // this should be something else, but now it is just
    // fraction good, i.e, linear decline
    return(fraction_good);
}

// fecundity function that accounts for costly resistance,
// as in Gandon & Vale eq (2)
double b(double const resistance)
{
    return(bmax * exp(-c * resistance));
}

// calculate the odds of a plasmid loss event
// and return this is in gamma val
// also determine whether this is a loss of a 
// good or a bad plasmid and return this in 
// loss_good
void gamma_loss(double &gamma_good, double &gamma_bad, Individual &host)
{
    assert(host.nplasmids_good <= host.nplasmids);

    gamma_good = host.nplasmids_good * gamma_G;
    gamma_bad = (host.nplasmids - host.nplasmids_good) * gamma_B;
}

// setup a distribution of events and choose
// which event to do
void event_chooser(int const time_step)
{
    // some auxiliary variables to store temporary values about rates
    double rate_birth, rate_loss_good, rate_loss_bad, rate_infect;

    // make a vector of total rates 
    // (these are the sums of all the individual 
    // infection, death and loss rates)
    //
    // we use this later to establish which of the events
    // will be picked
    int n_events = 10;
    std::vector <double> total_rates(n_events, 0.0);

    // declare a vector that contains all the rates
    // re an individual's birth
    std::vector <double> birth_rates_susceptible;
    std::vector <double> birth_rates_infected;

    // declare vectors containing infection rates of susceptibles
    // with good and bad plasmids respectively
    std::vector <double> infection_rate_susceptible_good;
    std::vector <double> infection_rate_susceptible_bad;
    
    // same for infected (aka rates as co-infection)
    std::vector <double> infection_rate_infected_good;
    std::vector <double> infection_rate_infected_bad;

    // declare vectors containing rates of
    // plasmid loss (good and bad respectively)
    std::vector <double> loss_rates_good;
    std::vector <double> loss_rates_bad;

    std::vector <double> death_rate_infected;
    
    double cumul_loss_rate_good = 0.0;
    double cumul_loss_rate_bad = 0.0;

    // go through all susceptibles and calculate birth
    // and infection rates
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        // 0. Birth of susceptibles
        rate_birth = b(Susceptible[S_idx].x) * (1.0 - kappa * N);
        birth_rates_susceptible.push_back(rate_birth);

        total_rates[0] += rate_birth;

        // 1. infection of susceptible by good plasmid
        // note that infection rates are independent of the number of 
        // bacteria already infected
        rate_infect = (1.0 - Susceptible[S_idx].x) * psi_G;
        infection_rate_susceptible_good.push_back(rate_infect);
        
        total_rates[1] += rate_infect;
        
        // 2. infection of susceptible by bad plasmid
        // note that infection rates are independent of the number of 
        // bacteria already infected
        rate_infect = (1.0 - Susceptible[S_idx].x) * psi_B;
        infection_rate_susceptible_bad.push_back(rate_infect);
        
        total_rates[2] += rate_infect;
    }

    // 3. Deaths susceptibles
    double death_rate_S = Ns * d;

    total_rates[3] += death_rate_S;

    double co_infection_bad, co_infection_good, death_rate;

    // now go through the infected hosts 
    for (int inf_idx = 0; inf_idx < Ni; ++inf_idx)
    {
        // 4. birth infected host
        rate_birth = F(Infected[inf_idx].fraction_good) * 
            b(Infected[inf_idx].x) * (1.0 - kappa * N);
        
        birth_rates_infected.push_back(rate_birth);

        total_rates[4] += rate_birth;

        // 5, 6. loss of a good or bad plasmid 
        // a single plasmid loss event
        gamma_loss(rate_loss_good, rate_loss_bad, Infected[inf_idx]);

        cumul_loss_rate_bad += rate_loss_bad;
        cumul_loss_rate_good += rate_loss_good;

        loss_rates_good.push_back(rate_loss_good);
        total_rates[5] += rate_loss_good;

        loss_rates_bad.push_back(rate_loss_bad);
        total_rates[6] += rate_loss_bad;

        // 7. death rate infected
        death_rate = Infected[inf_idx].fraction_good * dG + 
            (1.0 - Infected[inf_idx].fraction_good) * dB;

        death_rate_infected.push_back(death_rate);
        total_rates[7] += death_rate;

        // 9. infection with another plasmid
        co_infection_good = sigma * (1.0 - Infected[inf_idx].x) * psi_G;
        infection_rate_infected_good.push_back(co_infection_good);
        
        total_rates[8] += co_infection_good;

        co_infection_bad = sigma * (1.0 - Infected[inf_idx].x) * psi_B;
        infection_rate_infected_bad.push_back(co_infection_bad);
        
        total_rates[9] += co_infection_bad;
    } 

    // done, now determine what to do by making a weighted distribution
    // this will return a number between 0 and n_events - 1
    // dependent on the relative weighting of each event
    std::discrete_distribution<int> total_distribution(total_rates.begin(), total_rates.end());

    int S_idx,I_idx;

    int event_type = total_distribution(rng_r);

    // now sample a single sample from this discrete distribution
    // this is the type of event that will be chosen
    switch(event_type)
    {
        // next, we now determine the actual individuals
        // affected by the event
        
        case 0: // birth of susceptible
            {
            std::discrete_distribution <int> birth_susc_distribution(
                    birth_rates_susceptible.begin()
                    ,birth_rates_susceptible.end());

            S_idx = birth_susc_distribution(rng_r);

            assert(S_idx >= 0);
            assert(S_idx < Ns);

            // execute the birth() function which also updates the stats
            birth(Susceptible[S_idx], true);

            break;
            }

        case 1: // infection of a susceptible by good plasmid
            {
            std::discrete_distribution <int> infection_susceptible_good_dist(
                infection_rate_susceptible_good.begin()
                ,infection_rate_susceptible_good.end());

            S_idx = infection_susceptible_good_dist(rng_r);

            assert(S_idx >= 0);
            assert(S_idx < Ns);

            infection_susceptible(S_idx, true);

            break;
            }
        case 2: // infection of a susceptible by a bad plasmid
            {
            std::discrete_distribution <int> infection_susceptible_bad_dist(
                infection_rate_susceptible_bad.begin()
                ,infection_rate_susceptible_bad.end());

            S_idx = infection_susceptible_bad_dist(rng_r);

            assert(S_idx >= 0);
            assert(S_idx < Ns);

            infection_susceptible(S_idx, false);

            break;
            }
        case 3: // death susceptible
            death_susceptible();
            break;
        case 4: // birth infected host
            {
            std::discrete_distribution <int> birth_infected_dist(
                birth_rates_infected.begin()
                ,birth_rates_infected.end());

            I_idx = birth_infected_dist(rng_r);

            assert(I_idx >= 0);
            assert(I_idx < Ni);

            birth(Infected[I_idx], false);
            }
            break;

        case 5: // loss of a single good plasmid
            {
            std::discrete_distribution <int> loss_good_plasmid_dist(
                    loss_rates_good.begin()
                    ,loss_rates_good.end());

            I_idx = loss_good_plasmid_dist(rng_r);

            assert(I_idx >= 0);
            assert(I_idx < Ni);
            
            loss_plasmid(I_idx, true);
            
            break;
            }

        case 6: // loss of a single bad plasmid
            {
            std::discrete_distribution <int> loss_bad_plasmid_dist(
                    loss_rates_bad.begin()
                    ,loss_rates_bad.end());

            I_idx = loss_bad_plasmid_dist(rng_r);

            assert(I_idx >= 0);
            assert(I_idx < Ni);
            
            loss_plasmid(I_idx, false);
            
            break;
            }
        case 7:// death infected
            {
            std::discrete_distribution <int> death_infected_dist(
                    death_rate_infected.begin()
                    ,death_rate_infected.end());

            I_idx = death_infected_dist(rng_r);


            assert(I_idx >= 0);
            assert(I_idx < Ni);
            
            loss_plasmid(I_idx, false);

            death_infected(I_idx);
            }
        default:
            std::cout << "switch error" << std::endl;
            break;
    } // end switch
} // end event_chooser(int const time_step)

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
} // end main



