// sexual selection in bacteria

// initial model based on plasmid model by Gandon & Vale
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

int skip_output_rows = 10;

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
    // obtain all parameters from the command line
    max_time = atof(argv[1]);
    p_good_init = atof(argv[2]);
    kappa = atof(argv[3]);
    bmax = atof(argv[4]);
    c = atof(argv[5]);
    gamma_G = atof(argv[6]);
    gamma_B = atof(argv[7]);
    psi_G = atof(argv[8]);
    psi_B = atof(argv[9]);
    d = atof(argv[10]);
    dG = atof(argv[11]);
    dB = atof(argv[12]);
    sigma = atof(argv[13]);
    mu_x = atof(argv[14]);
    sdmu_x = atof(argv[15]);
    init_x = atof(argv[16]);
    base_name = argv[17];
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
        << "mu_x" << ";" << mu_x << std::endl
        << "sdmu_x" << ";" << sdmu_x << std::endl
        << "init_x" << ";" << init_x << std::endl
        << "n_plasmid_init" << ";" << n_plasmid_init << std::endl;

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

    assert(parent.x >= 0);
    assert(parent.x <= 1.0);

    // mutate the resistance allele
    kid.x = mutation(parent.x);

    // limit boundaries of resistance
    kid.x = std::clamp(kid.x, 0.0, 1.0);

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

    assert(random_susceptible >= 0);
    assert(random_susceptible < Ns);
    Susceptible[random_susceptible] = Susceptible[Ns - 1];
    --Ns;

    assert(Ns >= 0);
    assert(Ns + Ni <= N);
}// end death_susceptible()

// co infect a host (index I_idx)
// with another plasmid
void co_infection(int const I_idx, bool const plasmid_good)
{
    // check whether the individual is indeed infected
    assert(Infected[I_idx].nplasmids > 0);
    assert(Infected[I_idx].nplasmids <=nplasmid_max);
    // do nothing if individual already contains max number of
    // plasmids
    if (Infected[I_idx].nplasmids == nplasmid_max)
    {
        return;
    }

    // add another plasmid to the stack of plasmids
    Infected[I_idx].plasmid_good[Infected[I_idx].nplasmids++] = plasmid_good;

    // update the count of good plasmids
    // if the value of plasmid_good is false, it will add 0
    Infected[I_idx].nplasmids_good += plasmid_good;
    Infected[I_idx].fraction_good = (double)Infected[I_idx].nplasmids_good/Infected[I_idx].nplasmids;
} // end co_infection

// death of an infected individual at location I_idx
void death_infected(int const I_idx)
{
    assert(I_idx >= 0);
    assert(I_idx < Ni);

    assert(Infected[I_idx].nplasmids > 0);
    Infected[I_idx] = Infected[Ni - 1];
    --Ni;
}

void loss_plasmid(int const I_idx, bool const plasmid_good)
{
    assert(Infected[I_idx].nplasmids > 0);
    assert(I_idx >= 0);
    assert(I_idx < Ni);

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
    data_file << "time;mean_resistance;var_resistance;mean_plasmid_good;var_plasmid_good;Ns;Ni;mean_nplasmid;var_nplasmid;" << std::endl;
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

    if (Ns == 0 && Ni == 0)
    {
        exit(1);
    }

    // go through all susceptibles and calculate birth
    // and infection rates
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        
        if (Susceptible[S_idx].x < 0.0)
        {
            std::cout << "time: " << time_step << " Ns: " << Ns << " inf_idx: " << S_idx << " " << Susceptible[S_idx].x << std::endl;
        }

        assert(Susceptible[S_idx].x >= 0);
        assert(Susceptible[S_idx].x <= 1.0);

        // 0. Birth of susceptibles
        rate_birth = b(Susceptible[S_idx].x) * (1.0 - kappa * (Ns + Ni));
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
    } // end for S_idx

    // 3. Deaths susceptibles
    double death_rate_S = Ns * d;

    total_rates[3] += death_rate_S;

    double co_infection_bad, co_infection_good, death_rate;

    // now go through the infected hosts 
    for (int inf_idx = 0; inf_idx < Ni; ++inf_idx)
    {
        // bounds checking. If there is a buffer overflow
        // these numbers typically are not between bounds
        assert(Infected[inf_idx].x >= 0.0);
        assert(Infected[inf_idx].x <= 1.0);
        
        if (Infected[inf_idx].x < 0.0)
        {
            std::cout << "time: " << time_step << " Ni: " << Ni << " inf_idx: " << inf_idx << " " << Infected[inf_idx].x << std::endl;
        }


        // 4. birth infected host
        rate_birth = F(Infected[inf_idx].fraction_good) * 
            b(Infected[inf_idx].x) * (1.0 - kappa * (Ns + Ni));
        
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

        // 8. infection with a good plasmid
        co_infection_good = sigma * (1.0 - Infected[inf_idx].x) * psi_G;
        infection_rate_infected_good.push_back(co_infection_good);
        
        total_rates[8] += co_infection_good;

        // 9. infection with a bad plasmid
        co_infection_bad = sigma * (1.0 - Infected[inf_idx].x) * psi_B;
        infection_rate_infected_bad.push_back(co_infection_bad);
        
        total_rates[9] += co_infection_bad;
    }

    assert(infection_rate_infected_bad.size() == Ni);

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
        //
        // note additional curly braces in switch cases below,
        // coz https://stackoverflow.com/questions/92396/why-cant-variables-be-declared-in-a-switch-statement 
        
        case 0: // birth of susceptible
            {
                assert(birth_rates_susceptible.size() == Ns);
                // set up a probability distribution
                // that determines which individual will be drawn
                // to give birth
            std::discrete_distribution <int> birth_susc_distribution(
                    birth_rates_susceptible.begin()
                    ,birth_rates_susceptible.end());

            // then draw from the probability distribution
            // to determine the individual that gives birth
            S_idx = birth_susc_distribution(rng_r);

            // check whether the number is indeed
            // conforming to a susceptible
            //
            // note that S_idx is an index from 0 to Ns - 1
            // thus it does not include Ns itself
            assert(S_idx >= 0);
            assert(S_idx < Ns);

            // execute the birth() function which also updates the stats
            birth(Susceptible[S_idx], true);

            break;
            }

        case 1: // infection of a susceptible by good plasmid
            {
                assert(infection_rate_susceptible_good.size() == Ns);
                // set up a probability distribution
                // that determines which individual will get infected
                // by a good plasmid
                std::discrete_distribution <int> infection_susceptible_good_dist(
                    infection_rate_susceptible_good.begin()
                    ,infection_rate_susceptible_good.end());

                // then draw from the probability distribution
                // to determine the individual that gets infected
                S_idx = infection_susceptible_good_dist(rng_r);

                // bounds checking
                assert(S_idx >= 0);
                assert(S_idx < Ns);

                infection_susceptible(S_idx, true);

                break;
            }
        case 2: // infection of a susceptible by a bad plasmid
            {
                // set up a probability distribution
                // that determines which individual will get infected
                // by a bad plasmid
                std::discrete_distribution <int> infection_susceptible_bad_dist(
                    infection_rate_susceptible_bad.begin()
                    ,infection_rate_susceptible_bad.end());

                // then draw from the probability distribution
                // to determine the individual that gets infected
                S_idx = infection_susceptible_bad_dist(rng_r);

                assert(S_idx >= 0);
                assert(S_idx < Ns);

                infection_susceptible(S_idx, false);

                break;
            }
        case 3: // death susceptible
            // each susceptible has the same chance of dying
            // so we do not need to set up a probability distribution
            // determining which susceptible is more likely to die relative
            // to others, we simply pick a 
            // random individual

            // as we do not initialize any variables
            // within a 'case n:' clause
            // no curly braces needed
            death_susceptible();
            break;
        case 4: // birth infected host
            {
                // set up probability distribution that determines
                // which individual will give birth
                std::discrete_distribution <int> birth_infected_dist(
                    birth_rates_infected.begin()
                    ,birth_rates_infected.end());

                // draw an individual to give birth from that distribution
                I_idx = birth_infected_dist(rng_r);

                // bounds checking
                assert(I_idx >= 0);
                assert(I_idx < Ni);

                // perform a birth event
                birth(Infected[I_idx], false);
                
                break;
            }
        case 5: // loss of a single good plasmid
            {
                // set up probability distribution that determines
                // which individual will lose a good plasmid
                std::discrete_distribution <int> loss_good_plasmid_dist(
                        loss_rates_good.begin()
                        ,loss_rates_good.end());

                // draw an individual from that distribution which is going
                // to lose a plasmid
                I_idx = loss_good_plasmid_dist(rng_r);

                // bounds checking
                assert(I_idx >= 0);
                assert(I_idx < Ni);
               
                // perform the actual plasmid loss 
                loss_plasmid(I_idx, true);
            
                break;
            }
        case 6: // loss of a single bad plasmid
            {
                // set up probability distribution that determines
                // which individual will lose a bad plasmid
                std::discrete_distribution <int> loss_bad_plasmid_dist(
                        loss_rates_bad.begin()
                        ,loss_rates_bad.end());

                // draw an individual from that distribution which is going
                // to lose a plasmid
                I_idx = loss_bad_plasmid_dist(rng_r);

                // bounds checking
                assert(I_idx >= 0);
                assert(I_idx < Ni);

                // perform the actual plasmid loss
                loss_plasmid(I_idx, false);
                
                break;
            }
        case 7:// death infected
            {
                // set up probability distribution that determines
                // which individual will lose a bad plasmid
                std::discrete_distribution <int> death_infected_dist(
                        death_rate_infected.begin()
                        ,death_rate_infected.end());

                // draw the individual that is going to die from that distribution
                I_idx = death_infected_dist(rng_r);

                // bounds check
                assert(I_idx >= 0);
                assert(I_idx < Ni);
                
                // perform actual death 
                death_infected(I_idx);

                break;
            }
        case 8: // co-infection good plasmid
            {
                //set up probability distribution that determines
                //which individaul will get co-infected by a good plasmid
                std::discrete_distribution <int> co_infection_good_dist(
                    infection_rate_infected_good.begin()
                    ,infection_rate_infected_good.end());

                // draw the individual that is going to get co-infected
                I_idx = co_infection_good_dist(rng_r);

                // bounds check
                assert(I_idx >= 0);
                assert(I_idx < Ni);

                // perform the actual co-infection
                co_infection(I_idx, true);

                break;
            }
        case 9: // co-infection bad plasmid
            {
                //set up probability distribution that determines
                //which individaul will get co-infected by a good plasmid
                std::discrete_distribution <int> co_infection_bad_dist(
                    infection_rate_infected_bad.begin()
                    ,infection_rate_infected_bad.end());

                // draw the individual that is going to get co-infected
                I_idx = co_infection_bad_dist(rng_r);

                if (I_idx >= Ni)
                {
                    std::cout << I_idx << std::endl;
                }

                // bounds check
                assert(I_idx >= 0);
                assert(I_idx < Ni);

                // perform the actual co-infection
                co_infection(I_idx, false);

                break;
            }
        default:
            std::cout << "switch error" << std::endl;
            break;
    } // end switch
} // end event_chooser(int const time_step)

// write the data
void write_data(std::ofstream &data_file
        ,int const time_step)
{
    double mean_resistance = 0.0;
    double ss_resistance = 0.0;

    double mean_freq_plasmid_good = 0.0;
    double ss_freq_plasmid_good = 0.0;

    double mean_n_plasmid = 0.0;
    double ss_n_plasmid = 0.0;

    double x,p_good,n_plasmid;

    for (int I_idx = 0; I_idx < Ni; ++I_idx)
    {
        x = Infected[I_idx].x;

        mean_resistance += x;
        ss_resistance += x * x;

        p_good = Infected[I_idx].fraction_good;

        mean_freq_plasmid_good += p_good;
        ss_freq_plasmid_good += p_good * p_good;

        n_plasmid = Infected[I_idx].nplasmids;

        mean_n_plasmid += n_plasmid;
        ss_n_plasmid += n_plasmid * n_plasmid;
    }

    mean_resistance /= Ni;
    mean_freq_plasmid_good /= Ni;
    mean_n_plasmid /= Ni;

    double var_resistance = ss_resistance / Ni 
        - mean_resistance * mean_resistance;

    double var_freq_plasmid_good = ss_freq_plasmid_good / Ni
        - mean_freq_plasmid_good * mean_freq_plasmid_good;

    double var_n_plasmid = ss_n_plasmid / Ni
        - mean_n_plasmid * mean_n_plasmid;

    data_file << time_step << ";"
        << mean_resistance << ";"
        << var_resistance << ";"
        << mean_freq_plasmid_good << ";"
        << var_freq_plasmid_good << ";"
        << Ns << ";"
        << Ni << ";"
        << mean_n_plasmid << ";"
        << var_n_plasmid << ";" << std::endl;

}

int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    
    // initialize file to write data to
    std::string file_name = base_name + ".csv";
    std::ofstream data_file{file_name};

    // start with the parameters and then the data
    write_parameters(data_file);
    write_data_headers(data_file);

    // initialize the population
    init_pop();

    // numerically iterate the simulation
    for (int time_idx = 0; time_idx < max_time; ++time_idx)
    {
        event_chooser(time_idx);

        if (time_idx % skip_output_rows == 0)
        {
            write_data(data_file, time_idx);
        }
    }
    
} // end main



