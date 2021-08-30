// Fisherian sexual selection in bacteria

// initial model based on plasmid model by Gandon & Vale
// modified by Ana Duarte
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
#include "individual.hpp"


// C++ random number generation unsigned int seed = get_nanoseconds();
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r{seed};
std::uniform_real_distribution<> uniform(0.0,1.0);

// parameters are mostly place holders
// and will be set using init_arguments() on command line


// file base name
std::string base_name;

// Preference and trait alleles are boolean
// p == 0 is p1 (no preference), p == 1 is p2 (preference)
// t == 0 is t1 (no trait), t ==1 is t2 (trait present)

// initial frequency of preference allele
// initial frequency of trait allele
double init_p2 = 0;
double init_t2 = 0;

// frequency of p2 and t2 in the population
double freq_p2 = 0.0;
double freq_t2 = 0.0;

// max number of plasmids per cell
int const nplasmid_max = 1;

// number of timesteps that the simulation should run
int max_time = 10;

int skip_output_rows = 100;

// initial number of plasmids in infected cells
double n_plasmid_init = 1;

// initial frequency of uninfected (no plasmid)
double p_noplasmid_init = 1;

// density dependence
double kappa = 1;

// growth rate in absence of density dependence
double bmax = 1;

// preference for t == 1
double alpha = 0.0;

// preference cost factor 
// trait cost factor
double c = 0.0;
double epsilon = 0.0;

double cost_pref[3] = {0.0,0.0,0.0};
double cost_trait[3] = {0.0,0.0,0.0};

//cost of having plasmid
double delta = 0.0;

// loss rate of plasmid
double gamma_loss = 0.0;

// conjugation rate 
double pi = 0.0;

// vector to contain attractiveness
// of individuals to homozygous recipients
double attr_homz_recip[3] = {0.0,0.0,0.0};

// vector to contain attractiveness
// of individuals to heterozygous recipients
double attr_hetz_recip[3] = {0.0,0.0,0.0};
//
//
//chromosomal integration rate (plasmid integrated in chromosome)
// not in use at the moment
double tau = 0.0;

//plasmid formation rate (plasmid formed from chromosome)
// not in use currently
double lambda = 0.0;

// death rate 
double d = 0.0;

//recombination rate between plasmid and chromosome
double r = 0.0;

// preference, ornament mutation rate
// the rate at which p1 turns to p2 and vice-versa
// the rate at which t1 turns to t2 
double mu = 0.0;

// ornament biased mutation rate
double nu = 0.0;

// dominance coefficient 
// 0 = recessive, 0.5 = additive, 1 = dominant
double h = 0.0;  // for preference
double l = 0.0;  //for trait

// 4xS vectors of susceptible individuals (1 for each genotype)
std::vector < std::vector < Individual > > Susceptible;

// 4x4xI vectors of infected individuals 
std::vector < std::vector < std::vector < Individual > > > Infected;

// number of infected and susceptible hosts
//
int Ns, Ni, N;   

//// vector of attractiveness for all individuals 
//// from the perspective of a focal with homozygote preference 
//std::vector <double> attract_homozygote;
//
//// vector of attractiveness for all individuals 
//// from the perspective of a focal with a heterozygote preference 
//std::vector <double> attract_heterozygote;

// total birth rates of susceptible

// translate alleles to genotype
// 0: t1p1
// 1: t2p1
// 2: t1p2
// 3: t2p2
// first index: t
// second index: p
int alleles2genotypenr[2][2] = {{0,2},{1,3}};

bool geno_has_t2[4] = {0,1,0,1};
bool geno_has_p2[4] = {0,0,1,1};

int integer_division(int const x, int const y)
{
    return(int(floor(double(x) / y)));
}

// calculate attractiveness of Infected individual
// for a putative homozygote with preference
double calc_attract_susceptible(
        int const geno_susceptible
        ,int const geno_chr_infected
        ,int const geno_plasmid_infected)
{
    assert(geno_susceptible >= 0);
    assert(geno_susceptible < 4);
    assert(geno_chr_infected >= 0);
    assert(geno_chr_infected < 4);
    assert(geno_plasmid_infected >= 0);
    assert(geno_plasmid_infected < 4);

    // no preference, no need to calculate anything more
    if (!geno_has_p2[geno_susceptible])
    {
        return(1.0);
    }

    // check if infected has trait
	   // and susceptible has preference
    int n_trait_alleles = geno_has_t2[geno_chr_infected] + 
        geno_has_t2[geno_plasmid_infected]; 

    return(attr_homz_recip[n_trait_alleles]);

} // end calc_attractiveness

double calc_attract_infected(
        int const geno_chr_infected_recipient
        ,int const geno_plasmid_infected_recipient
        ,int const geno_chr_infected_donor
        ,int const geno_plasmid_infected_donor
        )
{
    assert(geno_chr_infected_donor >= 0);
    assert(geno_chr_infected_donor < 4);
    assert(geno_plasmid_infected_donor >= 0);
    assert(geno_plasmid_infected_donor < 4);
    assert(geno_chr_infected_recipient >= 0);
    assert(geno_chr_infected_recipient < 4);
    assert(geno_plasmid_infected_recipient >= 0);
    assert(geno_plasmid_infected_recipient < 4);

    int n_pref_alleles = geno_has_p2[geno_chr_infected_recipient] +
        geno_has_p2[geno_plasmid_infected_recipient];

    if (n_pref_alleles == 0)
    {
        return(1.0);
    }

    int n_trait_alleles = geno_has_t2[geno_chr_infected_donor] 
        + geno_has_t2[geno_plasmid_infected_donor]; 

    return(n_pref_alleles == 2 ?
            attr_homz_recip[n_trait_alleles]
            :
            attr_hetz_recip[n_trait_alleles]);
}

// initialize population
void init_pop()
{
    // calculate number of susceptibles
    Ns = round(N*p_noplasmid_init);  

    Ni = N - Ns;

    bool allele_chr_is_p2;
    bool allele_chr_is_t2;
    int genotype_chr; 

    // we need to initialize the 4 vectors for each
    // of the susceptible genotypes according to their 
    // frequencies (assuming no LD, that will be built up
    // stochastically anyway)

    std::vector <Individual> v;
    for (int i = 0; i < 4; ++i)
    {
        Susceptible.push_back(v);
    }

    for (int j = 0; j < 4; ++j)
    {
        // and initialize 4x4 infecteds
        // by re-using empty susceptible vector
        Infected.push_back(Susceptible);
    }

    // initialize susceptibles vectors
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        allele_chr_is_p2 = uniform(rng_r) < init_p2;
        allele_chr_is_t2 = uniform(rng_r) < init_t2;

        Individual init_ind;
        init_ind.p_chr = allele_chr_is_p2;
        init_ind.t_chr = allele_chr_is_t2;
        init_ind.has_plasmid = false;

        genotype_chr = alleles2genotypenr[allele_chr_is_t2][allele_chr_is_p2];

        Susceptible[genotype_chr].push_back(init_ind);
        assert(Susceptible[genotype_chr].size());
    }// for (int S_idx = 0; S_idx < Ns; ++S_idx)

    bool allele_plm_is_p2;
    bool allele_plm_is_t2;

    int genotype_plasmid;

    // initialize infected individuals
    for (int I_idx = 0; I_idx < Ni; ++I_idx)
    {
        allele_chr_is_p2 = uniform(rng_r) < init_p2;
        allele_chr_is_t2 = uniform(rng_r) < init_t2;
        allele_plm_is_p2 = uniform(rng_r) < init_p2;
        allele_plm_is_t2 = uniform(rng_r) < init_t2;

        Individual init_ind;
        init_ind.p_chr = allele_chr_is_p2;
        init_ind.t_chr = allele_chr_is_t2;
        init_ind.p_plasmid = allele_plm_is_p2;
        init_ind.t_plasmid = allele_plm_is_t2;
        
        init_ind.has_plasmid = true;

        genotype_chr = alleles2genotypenr[allele_chr_is_t2][allele_chr_is_p2];
        genotype_plasmid = alleles2genotypenr[allele_plm_is_t2][allele_plm_is_p2];

        Infected[genotype_chr][genotype_plasmid].push_back(init_ind);
    }

}//end void init_pop() 

void mutate_infected(int const genotype_chr
        ,int const genotype_plasmid)
{
    assert(genotype_chr >= 0);
    assert(genotype_chr < 4);
    assert(genotype_plasmid >= 0);
    assert(genotype_plasmid < 4);

    // mutation bias?
    double mu_t_chr = geno_has_t2[genotype_chr] ? nu : mu;
    double mu_t_plasmid = geno_has_t2[genotype_plasmid] ? nu : mu;

    double mu_p_ny[2] = {1.0 - mu, mu};
    double mu_t_chr_ny[2] = {1.0 - mu_t_chr, mu_t_chr};
    double mu_t_plasmid_ny[2] = {1.0 - mu_t_plasmid, mu_t_plasmid};

    double total_mu = 1.0 - (1.0 - mu)*(1.0 - mu) * (1.0 - mu_t_chr) * (1.0 - mu_t_plasmid);

    // draw a sample from cumulative distribution
    double random_cumul_sample = uniform(rng_r) * total_mu;

    double cumul_sum_mut = 0.0;

    bool done = false;

    // build cumulative distribution
    for (int mut_p_chr_idx = 0; mut_p_chr_idx < 2; ++mut_p_chr_idx)
    {
        for (int mut_t_chr_idx = 0; mut_t_chr_idx < 2; ++mut_t_chr_idx)
        {
            for (int mut_p_plasmid_idx = 0; mut_p_plasmid_idx < 2; ++mut_p_plasmid_idx)
            {
                for (int mut_t_plasmid_idx = 0; mut_t_plasmid_idx < 2; ++mut_t_plasmid_idx)
                {
                    if (mut_p_chr_idx == 0 && 
                            mut_t_chr_idx == 0 && 
                            mut_p_plasmid_idx == 0 && 
                            mut_t_plasmid_idx == 0)
                    {
                        continue;
                    }

                    cumul_sum_mut += mu_p_ny[mut_p_chr_idx] *
                        mu_t_chr_ny[mut_t_chr_idx] *
                        mu_p_ny[mut_p_plasmid_idx] *
                        mu_t_plasmid_ny[mut_t_plasmid_idx];

                    if (random_cumul_sample <= cumul_sum_mut && !done)
                    {
                        // actually mutate the thing
                        bool p_chr_old = geno_has_p2[genotype_chr];
                        bool p_plasmid_old = geno_has_p2[genotype_plasmid];

                        bool t_chr_old = geno_has_t2[genotype_chr];
                        bool t_plasmid_old = geno_has_t2[genotype_plasmid];

                        bool p_chr_new = mut_p_chr_idx == 0 ? p_chr_old : !p_chr_old;
                        bool p_plasmid_new = mut_p_plasmid_idx == 0 ? 
                            p_plasmid_old 
                            : 
                            !p_plasmid_old;

                        bool t_chr_new = mut_t_chr_idx == 0 ? t_chr_old : !t_chr_old;

                        bool t_plasmid_new = mut_t_plasmid_idx == 0 ? 
                            t_plasmid_old 
                            : 
                            !t_plasmid_old;

                        int new_genotype_chr = 
                            alleles2genotypenr[t_chr_new][p_chr_new];

                        int new_genotype_plasmid = 
                            alleles2genotypenr[t_plasmid_new][p_plasmid_new];

                        assert(Infected[genotype_chr][genotype_plasmid].size() > 0);
                        Infected[genotype_chr][genotype_plasmid].pop_back();

                        Individual ind;
                        ind.p_chr = p_chr_new;
                        ind.t_chr = t_chr_new;
                        ind.p_plasmid = p_plasmid_new;
                        ind.t_plasmid = t_plasmid_new;
                        ind.has_plasmid = true;
                        Infected[new_genotype_chr][new_genotype_plasmid].push_back(ind);

                        done = true;
                    }
                }
            }
        }
    }

    // total sum should be very similar
    assert(std::fabs(cumul_sum_mut - total_mu) < 1.0e-07);
} // end mutate_infected()

void mutate_susceptible(int const genotype)
{
    assert(genotype >= 0);
    assert(genotype < 4);

    std::vector<double> mutation_probs;

    // mutation bias?
    double mu_t = geno_has_t2[genotype] ? nu : mu;

    // option 1: mutate p but not t
    mutation_probs.push_back(mu * (1.0 - mu_t));
    // option 2: mutate t but not p
    mutation_probs.push_back((1.0 - mu) * mu_t);
    // option 3: mutate both
    mutation_probs.push_back(mu * mu_t);

    std::discrete_distribution <int> mutation_dist(
            mutation_probs.begin()
            ,mutation_probs.end());

    int mutation_type = mutation_dist(rng_r);

    assert(mutation_type >= 0);
    assert(mutation_type < 3);

    bool p_allele = geno_has_p2[genotype];
    bool t_allele = geno_has_t2[genotype];
    bool p_allele_new = p_allele;
    bool t_allele_new = t_allele;

    switch(mutation_type)
    {
        case 0: // mutate p, not t
            {
                p_allele_new = !p_allele;            
                t_allele_new = t_allele;            
                break;
            }
        case 1: // mutate t, not p
            {
                p_allele_new = p_allele;            
                t_allele_new = !t_allele;            
                break;
            }
        case 2:
            {
                p_allele_new = !p_allele;            
                t_allele_new = !t_allele;            
                break;
            }

        default:
            std::cout << "something went wrong..." << std::endl;
    }// end switch

    Individual mutated_ind = Susceptible[t_allele][p_allele];
    mutated_ind.p_chr = p_allele_new;
    mutated_ind.t_chr = t_allele_new;

    Susceptible[alleles2genotypenr[t_allele_new][p_allele_new]].push_back(mutated_ind);
                    
}

// genotype returns an int
// which can be used to update the genotype counts vectors
int genotype_haploid(bool const allele_p, bool const allele_t)
{
    if (allele_p == 0 && allele_t == 0) 
    {
        return(0);
    }

    if (allele_p == 0 && allele_t == 1) 
    {
        return(1);
    }
    
    if (allele_p == 1 && allele_t == 0) 
    {
        return(2);
    }

    return(3);
} // end of int genotype()

//int genotype_diploid(int const genotype_chr, int const genotype_pl)
//{
//     
//    if(genotype_chr == 0) {
//	return(genotype_pl);
//        }
//
//    if(genotype_chr == 1) {
//        return(genotype_pl + 4);
//        }
//    
//    if(genotype_chr == 2) {
//        return(genotype_pl + 8);
//        }
//
//    return(genotype_pl + 12);
//} // end of int genotype_diploid()

// initialize the parameters through the command line
void init_arguments(int argc, char ** argv)
{
    // obtain all parameters from the command line
    max_time = atof(argv[1]);
    p_noplasmid_init = atof(argv[2]);
    kappa = atof(argv[3]);
    bmax = atof(argv[4]);
    c = atof(argv[5]);
    epsilon = atof(argv[6]);
    delta = atof(argv[7]);
    gamma_loss = atof(argv[8]);
    pi = atof(argv[9]);
    tau = atof(argv[10]);
    lambda = atof(argv[11]);
    d = atof(argv[12]);
    r = atof(argv[13]);
    mu = atof(argv[14]);
    nu = atof(argv[15]);
    init_p2 = atof(argv[16]);
    init_t2 = atof(argv[17]);
    alpha = atof(argv[18]);
    h = atof(argv[19]);
    l = atof(argv[20]);
    base_name = argv[21];

    attr_homz_recip[0] = 1.0; // attractiveness individual without ornament
    attr_homz_recip[1] = 1.0 + l * alpha; // attractiveness individual with one ornament allele
    attr_homz_recip[2] = 1.0 + alpha; // attractiveness individual with one ornament allele
    
    attr_hetz_recip[0] = 1.0; // attractiveness individual without ornament
    attr_hetz_recip[1] = 1.0 + h * l * alpha; // attractiveness individual with one ornament allele
    attr_hetz_recip[2] = 1.0 + h * alpha; // attractiveness individual with one ornament allele

    cost_pref[0] = 0.0;
    cost_pref[1] = h * c;
    cost_pref[2] = c;

    cost_trait[0] = 0.0;
    cost_trait[1] = l * epsilon;
    cost_trait[2] = epsilon;

}//end init_arguments()

void write_parameters(std::ofstream &data_file)
{
    data_file << std::endl << std::endl
        << "N" << ";" << N << std::endl
        << "max_time" << ";" << max_time << std::endl
        << "p_noplasmid_init" << ";" << p_noplasmid_init << std::endl
        << "kappa" << ";" << kappa << std::endl
        << "bmax" << ";" << bmax << std::endl
        << "c" << ";" << c << std::endl
        << "alpha" << ";" << alpha << std::endl
	<< "epsilon" << ";" << epsilon << std::endl
        << "delta" << ";" << delta << std::endl
        << "gamma" << ";" << gamma_loss  << std::endl
        << "pi" << ";" << pi << std::endl
        << "tau" << ";" << tau << std::endl
        << "lambda" << ";" << lambda << std::endl
        << "d" << ";" << d << std::endl
        << "r" << ";" << r << std::endl
        << "seed" << ";" << seed << std::endl
        << "mu" << ";" << mu << std::endl
        << "nu" << ";" << nu << std::endl
        << "init_p2" << ";" << init_p2 << std::endl
        << "init_t2" << ";" << init_t2 << std::endl
        << "n_plasmid_init" << ";" << n_plasmid_init << std::endl
	<< "h" << ";" << h << std::endl
	<< "l" << ";" << l << std::endl;

} // end write_parameters()


// if mutation occurs 
// turn allele p1 into p2 and vice-versa
// same for trait locus
bool mutation(double const mu, bool const myallele)
{
    return(uniform(rng_r) < mu ? !myallele : myallele);
}

// infection event of a susceptible 
// by conjugation with infected individual 
void infection_susceptible(int const genotype_susceptible
        ,int const genotype_infected_chr
        ,int const genotype_infected_plasmid
        )
{
    assert(genotype_susceptible >= 0);
    assert(genotype_susceptible < 4);
    assert(genotype_infected_chr >= 0);
    assert(genotype_infected_chr < 4);
    assert(genotype_infected_plasmid >= 0);
    assert(genotype_infected_plasmid < 4);
    
    // and assign newly infected to stack of infected
    // individuals
    Individual new_infected;
    new_infected.p_chr = geno_has_p2[genotype_susceptible];
    new_infected.t_chr = geno_has_t2[genotype_susceptible];
    new_infected.p_plasmid = geno_has_p2[genotype_infected_plasmid];
    new_infected.t_plasmid = geno_has_t2[genotype_infected_plasmid];
    new_infected.has_plasmid = true;
    Infected[genotype_susceptible][genotype_infected_plasmid].push_back(new_infected);

    // remove random old susceptible 
    // by overwriting it with the susceptible
    // at the end of the stack of susceptibles
    Susceptible[genotype_susceptible].pop_back();
} //end infection_susceptible()

void conjugation_infected(
            int const recipient_chromosome
            ,int const recipient_plasmid
            ,int const donor_chromosome
            ,int const donor_plasmid
            )
{
    assert(recipient_chromosome>= 0);
    assert(recipient_plasmid>= 0);
    assert(donor_chromosome>= 0);
    assert(donor_plasmid>= 0);
    assert(recipient_chromosome < 4);
    assert(recipient_plasmid < 4);
    assert(donor_chromosome < 4);
    assert(donor_plasmid < 4);

    // infected host receives a plasmid next to one it already has
    // now we have to randomly choose between 2 plasmids
    // hence, only a change when 
    if (uniform(rng_r) < 0.5) 
    {
        Individual new_ind;
        new_ind.p_chr = geno_has_p2[recipient_chromosome];
        new_ind.t_chr = geno_has_t2[recipient_chromosome];
        new_ind.p_plasmid = geno_has_p2[donor_plasmid];
        new_ind.t_plasmid = geno_has_t2[donor_plasmid];
        new_ind.has_plasmid = true;

        Infected[recipient_chromosome][recipient_plasmid].pop_back();
        Infected[recipient_chromosome][donor_plasmid].push_back(new_ind);
	}
} // end of conjugation_infected()

// loss of a single plasmid
void loss_plasmid()
{
    // select random infected
    std::uniform_int_distribution <int> infected_dist(0, Ni - 1);

    int random_infected = infected_dist(rng_r);

    int sum_I = 0;

    for (int infected_chr_idx = 0;
            infected_chr_idx < 4; ++infected_chr_idx)
    {
        for (int infected_plasmid_idx = 0;
                infected_plasmid_idx < 4; ++infected_plasmid_idx)
        {
            sum_I += Infected[infected_chr_idx][infected_plasmid_idx].size();

            assert(Infected[infected_chr_idx][infected_plasmid_idx].size() <= Ni);

            if (random_infected <= sum_I)
            {
                assert(Infected[infected_chr_idx][infected_plasmid_idx].size() > 0);

                Individual new_ind(Susceptible[infected_chr_idx][0]);
                Susceptible[infected_chr_idx].push_back(new_ind);

                Infected[infected_chr_idx][infected_plasmid_idx].pop_back();
                return;
            }
        }
    }

    // this should never be reached
    std::cout << "should not happen: plasmid loss" << std::endl;
}// end loss_plasmid()

// recombination between plasmid and chromosome alleles
void recombination(Individual &ind)
{
    Individual old_ind = ind;

    if (uniform(rng_r) < r)
    {
        ind.p_plasmid = old_ind.p_chr;
        ind.p_chr = old_ind.p_plasmid;
        ind.t_plasmid = old_ind.t_chr;
        ind.t_chr = old_ind.t_plasmid;
    }
}// end void recombination()

// birth event of an individual
// is going to be different with 
void birth_susceptible(int const &genotype)
{
    Individual kid(Susceptible[genotype][0]);

    kid.has_plasmid = false;

    Susceptible[genotype].push_back(kid);
}

void birth_infected(int const &genotype_chr
        ,int const &genotype_plasmid)
{
    Individual kid(Infected[genotype_chr][genotype_plasmid][0]);

    assert(kid.has_plasmid);

    Infected[genotype_chr][genotype_plasmid].push_back(kid);
} // end of birth()

// death of susceptible individual
void death_susceptible()
{
    std::uniform_int_distribution <int> susceptible_dist(0, Ns - 1);
    int random_susceptible = susceptible_dist(rng_r);

    int sum_S = 0;

    for (int susceptible_chr_idx = 0;
            susceptible_chr_idx < 4; ++susceptible_chr_idx)
    {
        sum_S += Susceptible[susceptible_chr_idx].size();

        if (random_susceptible <= sum_S)
        {
            Susceptible[susceptible_chr_idx].pop_back();
            return;   
        }
    }

    std::cout << "this should not have happened: death susceptible" << std::endl;

    exit(1);
}// end death_susceptible()

// in this model the death rate for susceptible 
// and infected individuals is the same
void death_infected()
{
    std::uniform_int_distribution <int> infected_dist(0, Ni - 1);
    int random_infected = infected_dist(rng_r);

    int sum_I = 0;

    for (int infected_chr_idx = 0;
            infected_chr_idx < 4; ++infected_chr_idx)
    {
        for (int infected_plasmid_idx = 0;
                infected_plasmid_idx < 4; ++infected_plasmid_idx)
        {
            sum_I += Infected[infected_chr_idx][infected_plasmid_idx].size();

            assert(Infected[infected_chr_idx][infected_plasmid_idx].size() > 0);
            assert(Infected[infected_chr_idx][infected_plasmid_idx].size() <= Ni);

            if (random_infected <= sum_I)
            {
                assert(Infected[infected_chr_idx][infected_plasmid_idx].size() > 0);

                Infected[infected_chr_idx][infected_plasmid_idx].pop_back();
                return;
            }
        }
    }

    // this should never be reached
    std::cout << "should not happen: death_infected" << std::endl;
    exit(1);
} // end death_infected()

// death of an infected individual at location I_idx

// write headers to the datafile
void write_data_headers(std::ofstream &data_file)
{
    data_file << "time;";

    for (int genotype_chr_idx = 0; 
            genotype_chr_idx < 4; ++genotype_chr_idx) 
    {
        data_file 
            << "S" 
            << (geno_has_t2[genotype_chr_idx] ? "t2" : "t1") 
            << (geno_has_p2[genotype_chr_idx] ? "p2" : "p1") 
            << ";";

        for (int genotype_plasmid_idx = 0; 
                genotype_plasmid_idx < 4; ++genotype_plasmid_idx) 
        {
            data_file 
                << "I" 
                << (geno_has_t2[genotype_chr_idx] ? "t2" : "t1") 
                << (geno_has_p2[genotype_chr_idx] ? "p2" : "p1") 
                << (geno_has_t2[genotype_plasmid_idx] ? "t2" : "t1") 
                << (geno_has_p2[genotype_plasmid_idx] ? "p2" : "p1") 
                << ";";
        }
    }

    data_file
        << "mean_freq_p2_total;"
        << "mean_freq_p2_susceptible;"
        << "mean_freq_p2_infected;"
        << "mean_freq_p2_plasmid;"
        << "mean_freq_t2_total;"
        << "mean_freq_t2_susceptible;"
        << "mean_freq_t2_infected;"
        << "mean_freq_t2_plasmid;" << std::endl;
}// end write_data_headers()



// fecundity function that accounts for costly trait and preference,
// as in Gandon & Vale eq (2)
double b_Susceptible(bool const trait, bool const pref)
{
    return(bmax * exp(- c * pref - epsilon * trait));
}

// birth rate of an infected individual
double b_Infected(bool const trait_pl
        ,bool const pref_pl
        ,bool const trait_chr
        ,bool const pref_chr)
{
    int n_trait_alleles = trait_pl + trait_chr;
    int n_pref_alleles = pref_pl + trait_pl;

    return(bmax * exp(- cost_pref[n_pref_alleles] - cost_trait[n_trait_alleles] - delta));
} // end of b_Infected()

void update_counters()
{
    Ns = 0;
    Ni = 0;

    for (int geno_chr_idx = 0; geno_chr_idx < 4; ++geno_chr_idx)
    {
        Ns += Susceptible[geno_chr_idx].size();

        for (int geno_plasmid_idx = 0; geno_plasmid_idx < 4; ++geno_plasmid_idx)
        {
            Ni += Infected[geno_chr_idx][geno_plasmid_idx].size();
        }
    }

    N = Ns + Ni;
} // end update counters

// setup a distribution of events and choose
// which event to do
void event_chooser(int const time_step)
{
    update_counters();

    // make a vector of total rates 
    // (these are the sums of all the individual 
    // infection, death and loss rates)
    //
    // we use this later to establish which of the events
    // will be picked
    //
    // total number of events is: 
    // 0. birth susceptible, with 4 different birth rates, hence 4
    // 1. infection of 4 types of susceptible by 16 different types of bacteria
    //      resulting in 64 different types of infection events
    // 2. birth infected host (16 different types), 16
    // 3. loss of plasmid (1 event as long as gamma, 1
    // 4. conjugation between infected and infected (16 * 16 = 256)
    // 5. death susceptible, 1
    // 6. death infected, 1
    // 7. mutation susceptible
    // 8. mutation infected
    //      - in t2 at rate nu
    //      - in t1, p1 and p2 at rate mu
    //      - more likely in infecteds as they carry more loci 

    // we perform things into 2 stages: we first develop a cumulative
    // distribution split over 7 'main' events
    // once an event is chosen we can then select further into the number
    // of events fitting the categories
    int n_events = 9;

    // make a vector containing the total number of rates
    // from which we will eventually sample which event to choose
    // during each timestep
    std::vector <double> total_rates(n_events, 0.0);

    // declare a vector that contains all the rates
    // re birth of all genotypes. These are unidimensional
    // vectors, by using integer division and modulos we can
    // later recover the corresponding indices of the genotypes
    std::vector <double> birth_rates_susceptible;
    std::vector <double> birth_rates_infected;

    std::vector <double> mutation_susceptible;
    std::vector <double> mutation_infected;

    // make vector to store force of infection rates between
    // susceptible and infecteds
    std::vector <double> force_infection_susceptible;
    
    std::vector <double> force_infection_infected;
    
    if (Ns == 0 && Ni == 0)
    {
        std::cout << "time: " << time_step << "  --> Population extinct!" << std::endl;
        exit(1);
    }

    if (Ni == 0)
    {
        std::cout << "time: " << time_step << "  --> Ni extinct" << std::endl;
        exit(1);
    }

    // aux variable to assess whether susceptible
    // recipient is choosy
    bool recipient_choosy;

    // and in case recipient is infected how many choosiness alleles
    // it actually has
    int n_choosiness_alleles;

    // aux variable to have mass-action in the infection rates
    // this is according to frequency-dependence/mass action 
    // (p17 in Keeling & Rohani)
    double inv_popsize = 1.0 / N;

    // aux variable for the force of infection
    double force_infection_ij;

    // density dependent correction of growth rates
    double dens_dep = 1.0 - kappa * (Ns + Ni);

    bool allele_chr_is_t2;
    bool allele_chr_is_p2;

    double birth_rate_i;

    double total_mu;

    // go through all susceptibles and calculate birth
    // and infection rates
    for (int geno_sus_chr_idx = 0; geno_sus_chr_idx < 4; ++geno_sus_chr_idx)
    {
        // 0. birth rates infected
        birth_rate_i = dens_dep <= 0.0 ?
            0.0
            :
            b_Susceptible(
                    geno_has_t2[geno_sus_chr_idx]
                    ,geno_has_p2[geno_sus_chr_idx]) * dens_dep * 
                        Susceptible[geno_sus_chr_idx].size();

        // add this birth rate to the stack
        birth_rates_susceptible.push_back(birth_rate_i);

        total_rates[0] += birth_rate_i;


        // 1. infection of susceptible with plasmid by 
        // infected individual
        //
        // iterate over all infected genotypes (chromosome x plasmid)
        for (int geno_inf_chr_idx = 0; 
                geno_inf_chr_idx < 4; ++geno_inf_chr_idx)
        {
            for (int geno_inf_plasmid_idx = 0; 
                    geno_inf_plasmid_idx < 4; ++geno_inf_plasmid_idx)
            {
                force_infection_ij = (1.0 - pi) * 
                           calc_attract_susceptible(
                                geno_sus_chr_idx
                               ,geno_inf_chr_idx 
                               ,geno_inf_plasmid_idx) 
                           * inv_popsize * Susceptible[geno_sus_chr_idx].size() * 
                                Infected[geno_inf_chr_idx][geno_inf_plasmid_idx].size();
            
                force_infection_susceptible.push_back(force_infection_ij);
                
                total_rates[1] += force_infection_ij;
            }
        }
        
        // 7. mutation susceptible
        // total mutation rate of this genotype
        // p always normal mutation rate
        // t2 may have higher mutation rate
        //
        //
        // probability of at least 1 mutation is 
        // 1 - (1-mu_p)*(1-mu_t)
        total_mu = 1.0 - (1.0 - mu) * 
            (1.0 - (geno_has_t2[geno_sus_chr_idx] ? nu : mu)) * 
                Susceptible[geno_sus_chr_idx].size();
        
        mutation_susceptible.push_back(total_mu);

        total_rates[7] += total_mu; 
    } // end for  for (int geno_sus_chr_idx = 0; geno_sus_chr_idx < 4; ++geno_sus_chr_idx)
    
    int geno_ctr = 0;

    // now go through the infected host genotypes
    for (int geno_inf_chr_idx = 0; 
            geno_inf_chr_idx < 4; ++geno_inf_chr_idx)
    {
        for (int geno_inf_plasmid_idx = 0; 
                geno_inf_plasmid_idx < 4; ++geno_inf_plasmid_idx)
        {
            // 2. birth infected host   
            birth_rate_i = dens_dep <= 0.0 ? 
                0.0 : 
                b_Infected(
                    geno_has_t2[geno_inf_plasmid_idx]
                    ,geno_has_p2[geno_inf_plasmid_idx]
                    ,geno_has_t2[geno_inf_chr_idx]
                    ,geno_has_p2[geno_inf_chr_idx]
                    ) * dens_dep * 
                        Infected[geno_inf_chr_idx][geno_inf_chr_idx].size();
        
            birth_rates_infected.push_back(birth_rate_i);

            total_rates[2] += birth_rate_i;
            
            // 8. mutation infected
            // total mutation rate of this genotype
            // p always normal mutation rate
            // t2 may have higher mutation rate
            //
            // probability of at least 1 mutation 
            // 1 - (1-mu_p)^2*(1-mu_t)^2
            total_mu = 1.0 - 
                    (1.0 - mu)*(1.0-mu) * 
                        (1.0 - (geno_has_t2[geno_inf_chr_idx] ? nu : mu)) *
                        (1.0 - (geno_has_t2[geno_inf_plasmid_idx] ? nu : mu)) *
                            Infected[geno_inf_chr_idx][geno_inf_plasmid_idx].size();

            mutation_infected.push_back(total_mu);

            total_rates[8] += total_mu;
        } // geno_inf_plasmid_idx
    } // end for geno_inf_chr_idx

    // 3. loss of plasmid
    total_rates[3] += gamma_loss * Ns;

    // 4. conjugation between infected and infected
    // 16 x 16 = 256 combinations
    for (int geno_inf_recip_chr_idx = 0; 
            geno_inf_recip_chr_idx < 4; ++geno_inf_recip_chr_idx)
    {
        for (int geno_inf_recip_plasmid_idx = 0; 
                geno_inf_recip_plasmid_idx < 4; ++geno_inf_recip_plasmid_idx)
        {
            for (int geno_inf_donor_chr_idx = 0; 
                    geno_inf_donor_chr_idx < 4; ++geno_inf_donor_chr_idx)
            {
                for (int geno_inf_donor_plasmid_idx = 0; 
                        geno_inf_donor_plasmid_idx < 4; ++geno_inf_donor_plasmid_idx)
                {
                    force_infection_ij = (1.0 - pi) * inv_popsize *
                        calc_attract_infected(geno_inf_recip_chr_idx
                                ,geno_inf_recip_plasmid_idx
                                ,geno_inf_donor_chr_idx
                                ,geno_inf_donor_plasmid_idx) * 
                                    Infected[geno_inf_recip_plasmid_idx][geno_inf_recip_plasmid_idx].size() *

                                    Infected[geno_inf_donor_plasmid_idx][geno_inf_donor_plasmid_idx].size();

                    // add infection force for this pair to the stack
                    force_infection_infected.push_back(force_infection_ij);
                    total_rates[4] += force_infection_ij;
                }
            }
        }
    } // end for int inf_idx;

    // 5. Deaths susceptibles
    double death_rate_S = Ns * d;

    total_rates[5] = death_rate_S;

    // 6. death rate infected
    double death_rate_I = Ni * d; 

    total_rates[6] = death_rate_I;

    // done, now determine what to do by making a weighted distribution
    // this will return a number between 0 and n_events - 1
    // dependent on the relative weighting of each event
    std::discrete_distribution<int> total_distribution(total_rates.begin(), total_rates.end());

    // sample from distribution
    int event_type = total_distribution(rng_r);

    int susceptible_genotype_idx;

    int SI_genotype_pair_idx;
    int genotype_infected_plasmid, genotype_infected_chr, genotype_susceptible;
    // now sample a single sample from this discrete distribution
    // this is the type of event that will be chosen
    switch(event_type)
    {
        // next, we now determine the actual individuals
        // affected by the event
        
        case 0: // birth of susceptible
        {
            // set up a probability distribution
            // that determines which individual will be drawn
            // to give birth
            std::discrete_distribution <int> birth_susc_distribution(
                birth_rates_susceptible.begin()
                ,birth_rates_susceptible.end());

            // then draw from the probability distribution
            // to determine the individual that gives birth
            susceptible_genotype_idx = birth_susc_distribution(rng_r);

            // note that S_idx is an index from 0 to Ns - 1
            // thus it does not include Ns itself
            assert(susceptible_genotype_idx >= 0);
            assert(susceptible_genotype_idx < 4);

            // execute the birth() function which also updates the stats
            birth_susceptible(susceptible_genotype_idx);

            break;
        } // end case 0

        case 1: // infection of a susceptible 
        {
            // choose a pair of susceptible and infected
            // individuals and tranmsmit plasmid
            //
            // first generate the sampling distribution from the infection rates
            std::discrete_distribution <int> infection_susceptible_dist(
                force_infection_susceptible.begin()
                ,force_infection_susceptible.end()
                );

            // obtain the pair of infected and susceptible
            SI_genotype_pair_idx = infection_susceptible_dist(rng_r);

            // get infected plasmid 
            genotype_infected_plasmid = SI_genotype_pair_idx % 4;

            // get infected chromosome 
            genotype_infected_chr = integer_division(SI_genotype_pair_idx,4) % 4;

            // get susceptible chromosome 
            genotype_susceptible = integer_division(
                    integer_division(
                        SI_genotype_pair_idx,4),4) % 4;

            // obtain susceptible idx
            infection_susceptible(
                    genotype_susceptible
                    ,genotype_infected_chr
                    ,genotype_infected_plasmid
                    );
            break;
        }

        case 2: // birth infected host
        {
            // set up probability distribution that determines
            // which individual will give birth
            std::discrete_distribution <int> birth_infected_dist(
                birth_rates_infected.begin()
                ,birth_rates_infected.end());

            // draw an individual to give birth from that distribution
            int birth_infected_idx = birth_infected_dist(rng_r);

            // this is a number that reflects a 2-dimensional array
            // so we need to disentangle which plasmid genotype it has
            genotype_infected_plasmid = birth_infected_idx % 4;

            // and which chromosomal genotype it has
            genotype_infected_chr = integer_division(birth_infected_idx,4) % 4;

            // perform a birth event
            birth_infected(
                    genotype_infected_chr
                    ,genotype_infected_plasmid);
            
            break;
        } // end case 2

        case 3: // loss of a single plasmid
        {
            // perform the actual plasmid loss 
            loss_plasmid();
        
            break;
        }

        case 4: // conjugation infected 
        {
            //randomly pick
            //which individual will be a receiver 
            // draw the individual that is going to get co-infected
            std::discrete_distribution <int> infection_infected_dist(
                force_infection_infected.begin()
                ,force_infection_infected.end());

            int II_pair_idx = infection_infected_dist(rng_r);

            // get genotype of donor plasmid
            int donor_plasmid = II_pair_idx % 4;

            int donor_chromosome = integer_division(II_pair_idx,4) % 4;

            int recipient_plasmid = 
                integer_division(
                    integer_division(II_pair_idx, 4)
                    ,4) % 4;

            int recipient_chromosome = 
                integer_division(
                    integer_division(
                        integer_division(
                            II_pair_idx, 4)
                            ,4)
                        ,4
                    ) % 4;

            // perform the actual conjugation
            // where a donor will be selected 
            // according to the receiver's preference
            conjugation_infected(
                    recipient_chromosome
                    ,recipient_plasmid
                    ,donor_chromosome
                    ,donor_plasmid
                    );
            break;
        }

        case 5: // death susceptible
	    {
            death_susceptible();

            break;
	    } // end case 5

        case 6:// death infected
        {
            // perform actual death 
            death_infected();

            break;
        }
        
        // mutational event in an infected individual 
        // (more likely coz 2 genomes)
        case 7:
        {
            std::discrete_distribution <int> mutation_susceptible_dist(
                mutation_susceptible.begin()
                ,mutation_susceptible.end());

            susceptible_genotype_idx = mutation_susceptible_dist(rng_r);

            mutate_susceptible(susceptible_genotype_idx);
        }

        case 8:
        {
            std::discrete_distribution <int> mutation_infected_dist(
                mutation_infected.begin()
                ,mutation_infected.end());

            int infected_genotype_idx = mutation_infected_dist(rng_r);

            int infected_genotype_plasmid = infected_genotype_idx % 4;
            int infected_genotype_chr = integer_division(infected_genotype_idx, 4) % 4;

            mutate_infected(infected_genotype_chr, infected_genotype_plasmid);
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
    // calculate allele freqs
    double mean_freq_p2_total = 0;
    double mean_freq_p2_susceptible = 0;
    double mean_freq_p2_infected = 0;
    double mean_freq_p2_chr = 0;
    double mean_freq_p2_plasmid = 0;

    double mean_freq_t2_total = 0;
    double mean_freq_t2_susceptible = 0;
    double mean_freq_t2_infected = 0;
    double mean_freq_t2_chr = 0;
    double mean_freq_t2_plasmid = 0;

    //double p2, p2_chr, p2_plasmid, t2, t2_chr, t2_plasmid, n_plasmid;

    data_file << time_step << ";";

    for (int genotype_chr_idx = 0; 
            genotype_chr_idx < 4; ++genotype_chr_idx) 
    {
        data_file << Susceptible[genotype_chr_idx].size() << ";";

        // calculate allele freqs
        mean_freq_p2_total += geno_has_p2[genotype_chr_idx] * 
            Susceptible[genotype_chr_idx].size();

        mean_freq_p2_susceptible += geno_has_p2[genotype_chr_idx] *
            Susceptible[genotype_chr_idx].size();

        mean_freq_p2_chr += geno_has_p2[genotype_chr_idx] *
            Susceptible[genotype_chr_idx].size();
       

        mean_freq_t2_total += geno_has_t2[genotype_chr_idx] *
            Susceptible[genotype_chr_idx].size();

        mean_freq_t2_susceptible += geno_has_t2[genotype_chr_idx] *
            Susceptible[genotype_chr_idx].size();
        
        mean_freq_t2_chr += geno_has_t2[genotype_chr_idx] *
            Susceptible[genotype_chr_idx].size();


        for (int genotype_plasmid_idx = 0;
                genotype_plasmid_idx < 4; ++genotype_plasmid_idx)
        {
            data_file << Infected[genotype_chr_idx][genotype_plasmid_idx].size() << ";";
            
            // update allele freqs for the chromosome of infected: p
            mean_freq_p2_total += geno_has_p2[genotype_chr_idx] * 
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_p2_chr += geno_has_p2[genotype_chr_idx] * 
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_p2_infected += geno_has_p2[genotype_chr_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();
            
            // update allele freqs for the plasmid of infected: p
            mean_freq_p2_total += geno_has_p2[genotype_plasmid_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_p2_plasmid += geno_has_p2[genotype_plasmid_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_p2_infected += geno_has_p2[genotype_plasmid_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            // update allele freqs for the chromosome of infected: t
            mean_freq_t2_total += geno_has_t2[genotype_chr_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_t2_chr += geno_has_t2[genotype_chr_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_t2_infected += geno_has_t2[genotype_chr_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            // update allele freqs for the plasmid of infected: t
            mean_freq_t2_total += geno_has_t2[genotype_plasmid_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_t2_plasmid += geno_has_t2[genotype_plasmid_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();

            mean_freq_t2_infected += geno_has_t2[genotype_plasmid_idx] *
                Infected[genotype_chr_idx][genotype_plasmid_idx].size();
        }
    }

    assert(mean_freq_p2_total <= Ns + 2 * Ni);
    assert(mean_freq_p2_susceptible <= Ns);
    assert(mean_freq_p2_infected <= 2*Ni);
    assert(mean_freq_p2_chr <= Ni + Ns);
    assert(mean_freq_p2_plasmid <= Ni);

    assert(mean_freq_t2_total <= Ns + 2 * Ni);
    assert(mean_freq_t2_susceptible <= Ns);
    assert(mean_freq_t2_infected <= 2*Ni);
    assert(mean_freq_t2_chr <= Ni + Ns);
    assert(mean_freq_t2_plasmid <= Ni);
    
    mean_freq_p2_total /= Ns + 2 * Ni;
    mean_freq_p2_susceptible /= Ns;
    mean_freq_p2_infected /= 2*Ni;
    mean_freq_p2_plasmid /= Ni;
    mean_freq_p2_chr /= Ni + Ns;

    mean_freq_t2_total /= Ns + 2 * Ni;
    mean_freq_t2_susceptible /= Ns;
    mean_freq_t2_infected /= 2*Ni;
    mean_freq_t2_plasmid /= Ni;
    mean_freq_t2_chr /= Ni + Ns;

	data_file 
        << mean_freq_p2_total << ";"
        << mean_freq_p2_susceptible << ";"
        << mean_freq_p2_infected << ";"
        << mean_freq_p2_plasmid << ";"
        << mean_freq_p2_chr << ";"
        << mean_freq_t2_total << ";"
        << mean_freq_t2_susceptible << ";"
        << mean_freq_t2_infected << ";"
        << mean_freq_t2_plasmid << ";" 
        << mean_freq_t2_chr << ";" 
        << std::endl;
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
        // I would not reset every timestep, rather
        // just update counts (easier)
        //reset_genotype_counts();
        std::cout << time_idx << std::endl;

        event_chooser(time_idx);

        if (time_idx % skip_output_rows == 0)
        {
            write_data(data_file, time_idx);
        }
    }
    
} // end main



