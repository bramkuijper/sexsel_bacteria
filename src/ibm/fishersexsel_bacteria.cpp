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


// C++ random number generation unsigned int seed = get_nanoseconds();
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r{seed};
std::uniform_real_distribution<> uniform(0.0,1.0);

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

// population size
int const N = 10000;

int const nplasmid_max = 1;

// number of timesteps that the simulation should run
int max_time = 10;

int skip_output_rows = 10;

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

//cost of having plasmid
double delta = 0.0;

// loss rate of plasmid
double loss_gamma = 0.0;

// conjugation rate 
double psi = 0.0;

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

// preference mutation rate
// the rate at which p1 turns to p2 and vice-versa
double mu_p = 0.0;

// trait mutation rate
// the rate at which t1 turns to t2 and vice-versa
double mu_t = 0.0;

// dominance coefficient 
// 0 = recessive, 0.5 = additive, 1 = dominant
double h = 0.0;  // for preference
double l = 0.0;  //for trait

// individual
struct Individual
{
    // genotype	
    //preference in chromosome 
    //trait in chromosome
    bool p_chr;
    bool t_chr;

    //preference in plasmid 
    //trait in plasmid
    bool p_plasmid;
    bool t_plasmid;

    // phenotype  ===> do I really need this? 
    double p_phen;
    double t_phen;
    
    // list of plasmids
    // note that we allow potentially for 
    // larger number of plasmids / cell than just 1 
    bool plasmid[nplasmid_max];

    // actual number of plasmids carried by individual
    int nplasmids;

    // count of the number of  plasmids
    // needs to be recalculated once we make changes to
    // plasmid composition
    int nplasmids_count;

};

typedef Individual Population[N];

Population Susceptible;
Population Infected;

// number of infected and susceptible hosts
//
int Ns, Ni;   

// vector of attractiveness for individuals with homozygote preference 
// vector of attractiveness for individuals with heterozygote preference 
std::vector <double> attract_homozygote;
std::vector <double> attract_heterozygote;
// for susceptible and infected

// calculate attractiveness of Infected individual
// for a putative homozygote with preference
double calc_attract_homozygote(int const Donor_idx)
{
    assert(Infected[Donor_idx].nplasmids > 0);
    assert(Infected[Donor_idx].nplasmids <= nplasmid_max);
    // check if infected has trait
	   // and susceptible has preference
    bool hastrait = Infected[Donor_idx].t_plasmid || Infected[Donor_idx].t_chr;
    double l_current; // trait dominance coeff 

    // if infected has trait
    // calculate its attractiveness
    if(hastrait)
    {	
	    // if heterozygous, use parameter dominance coefficient
	if(Infected[Donor_idx].t_plasmid != Infected[Donor_idx].t_chr)
	{
	    l_current = l;
	} else {
		   l_current = 1;
	       }
	// return attractiveness - alpha times the dominance coeff
	return(1.0 + l_current * alpha);

    } else { 
	       return(1.0);
           }
} // end calc_attractiveness

// calculate attractiveness of Infected with trait  
// for a putative heterozygote with preference (must be Infected if hetero)  
double calc_attract_heterozygote(int const Donor_idx)
{
    assert(Infected[Donor_idx].nplasmids > 0);
    assert(Infected[Donor_idx].nplasmids <= nplasmid_max);
    // check if infected has trait
    bool hastrait = Infected[Donor_idx].t_plasmid || Infected[Donor_idx].t_chr;
    double l_current; // trait dominance coeff 

    if(hastrait)
    {	
	    // if donor heterozygous, use parameter dominance coefficient
	if(Infected[Donor_idx].t_plasmid != Infected[Donor_idx].t_chr)
	{
	    l_current = l;
	} else { 
	           l_current = 1;
	       }
	// return attractiveness - alpha times the dominance coeff
	return(1.0 + h * l_current * alpha);

    } else { 
	       return(1.0);
    	   }
	  
} // end calc_attractiveness

// initialize population
void init_pop()
{
    
    Ns = round(N*p_noplasmid_init);  

    Ni = N - Ns;

    std::cout << "Initializing population...  Ns: " << Ns << " Ni: " << Ni << " " << std::endl;
    // initialize susceptibles
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
        Susceptible[S_idx].p_chr = uniform(rng_r) < init_p2;
        Susceptible[S_idx].t_chr = uniform(rng_r) < init_t2;
        Susceptible[S_idx].nplasmids = 0;
        Susceptible[S_idx].plasmid[0] = false;
    }// for (int S_idx = 0; S_idx < Ns; ++S_idx)

    // initialize infected individuals
    for (int I_idx = 0; I_idx < Ni; ++I_idx)
    {
        Infected[I_idx].p_chr = uniform(rng_r) < init_p2;
        Infected[I_idx].t_chr = uniform(rng_r) < init_t2;
        Infected[I_idx].p_plasmid = uniform(rng_r) < init_p2;
        Infected[I_idx].t_plasmid = uniform(rng_r) < init_t2;
        Infected[I_idx].nplasmids = n_plasmid_init;
        Infected[I_idx].plasmid[0] = true;
	attract_homozygote.push_back(calc_attract_homozygote(I_idx));
	attract_heterozygote.push_back(calc_attract_heterozygote(I_idx));

    }
   
    assert(attract_homozygote.size() == Ni); 
    assert(attract_heterozygote.size() == Ni); 
    
}//end void init_pop() 

// updating attractiveness vectors for birth/death events

// update attractiveness vectors when individual is dead
// only Infected individuals, Susceptible do not act as donors 
// hence have no attractiveness
void update_death_attract(int const dead_idx)  
{
    // update in the same way as dead individuals are updated in population
    // copy last individual on the stack over dead individual
    // then remove last element of vector
// this needs to happen before actual removal/addition of individuals
// to population stack
// hence, when asserting vector sizes we must account for the individual we just added/removed
// as population sizes will be updated after updating attractiveness 
    attract_homozygote[dead_idx] = attract_homozygote.back();
    attract_homozygote.pop_back();	   
    assert(attract_homozygote.size() == Ni - 1);

    attract_heterozygote[dead_idx] = attract_heterozygote.back();
    attract_heterozygote.pop_back();	   
    assert(attract_heterozygote.size() == Ni - 1);
}

// update attractiveness vectors when 
// infected individual reproduces 
// or a susceptible is newly infected  
void update_birth_attract(int const new_idx)  
{
    // birth/creation of an infected individual
    // update in the same way as births are updated in population
    // add new individual to end of attractiveness vector 
    attract_homozygote.push_back(calc_attract_homozygote(new_idx));
    attract_heterozygote.push_back(calc_attract_heterozygote(new_idx));
    assert(attract_homozygote.size() == Ni);
    assert(attract_heterozygote.size() == Ni);
}

    // conjugation between two infected individuals
void update_conj_attract(int const receiver_idx)  
{
    // recalculate individual's attraction value 
    // for putative homozygote and putative heterozygote 
    attract_homozygote[receiver_idx] = calc_attract_homozygote(receiver_idx);
    attract_heterozygote[receiver_idx] = calc_attract_heterozygote(receiver_idx);

    assert(attract_homozygote.size() == Ni);
    assert(attract_heterozygote.size() == Ni); 
     
} // end of update_conj_attract

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
    loss_gamma = atof(argv[8]);
    psi = atof(argv[9]);
    tau = atof(argv[10]);
    lambda = atof(argv[11]);
    d = atof(argv[12]);
    r = atof(argv[13]);
    mu_p = atof(argv[14]);
    mu_t = atof(argv[15]);
    init_p2 = atof(argv[16]);
    init_t2 = atof(argv[17]);
    alpha = atof(argv[18]);
    h = atof(argv[19]);
    l = atof(argv[20]);
    base_name = argv[21];
    
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
        << "loss_gamma" << ";" << loss_gamma << std::endl
        << "psi" << ";" << psi << std::endl
        << "tau" << ";" << tau << std::endl
        << "lambda" << ";" << lambda << std::endl
        << "d" << ";" << d << std::endl
        << "r" << ";" << r << std::endl
        << "seed" << ";" << seed << std::endl
        << "mu_p" << ";" << mu_p << std::endl
        << "mu_t" << ";" << mu_t << std::endl
        << "init_p2" << ";" << init_p2 << std::endl
        << "init_t2" << ";" << init_t2 << std::endl
        << "n_plasmid_init" << ";" << n_plasmid_init << std::endl
	<< "h" << ";" << h << std::endl
	<< "l" << ";" << l << std::endl;

} // end write_parameters()


// if mutation occurs 
// turn allele p1 into p2 and vice-versa
// same for trait locus
double mutation(double const mu, bool const myallele)
{
    bool allele = myallele; 
    if (uniform(rng_r) < mu)
    {
	allele = !allele;
        return(allele);
    }

    return(allele);
}

// infection event of a susceptible 
// by conjugation with infected individual 
void infection_susceptible(int const S_idx)
{
    assert(S_idx >= 0);
    assert(S_idx < Ns);
    
    // check indeed that the susceptible individual
    // is not infected
    assert(Susceptible[S_idx].nplasmids == 0);

    int I_idx;
    // CHOOSING AN INFECTED INDIVIDUAL
    // if susceptible has no preference
    // choose an infected individual randomly 
    if(!Susceptible[S_idx].p_chr)  
    {
        std::uniform_int_distribution<int> infected_sampler(0, Ni - 1);
	I_idx = infected_sampler(rng_r);
    }

    // if susceptible has preference
    // generate weighted distribution with attractiveness values
    // of infected individuals
    // for a putative homozygote
    if(Susceptible[S_idx].p_chr)
    {
        std::discrete_distribution <int> infected_attract_dist(attract_homozygote.begin(), attract_homozygote.end());
	I_idx = infected_attract_dist(rng_r);
    }

    assert(I_idx >= 0);
    assert(I_idx < Ni);

    // assert the infected IS infected

    assert(Infected[I_idx].nplasmids > 0);
    assert(Infected[I_idx].nplasmids <= nplasmid_max);

    // and assign newly infected to stack of infected
    // individuals

    Infected[Ni] = Susceptible[S_idx];

    // remove old susceptible 
    // by overwriting it with the susceptible
    // at the end of the stack of susceptibles
    Susceptible[S_idx] = Susceptible[Ns - 1];
    --Ns;

    //copy plasmid genotype of donor to receiver cell
    Infected[Ni].p_plasmid = Infected[I_idx].p_plasmid;
    Infected[Ni].t_plasmid = Infected[I_idx].t_plasmid;

    //update boolean and count of plasmids -- keeping this in case it's useful 
    Infected[Ni].plasmid[0] = 1;  // this works because the newly infected cell had no plasmids before
    ++Infected[Ni].nplasmids; 

    //boundary checks
    assert(Infected[Ni].nplasmids > 0);
    assert(Infected[Ni].nplasmids <= nplasmid_max);
     
    ++Ni;
    //update attractiveness vectors
    // there's a new infected in town!
    update_birth_attract(Ni - 1);

    // check whether population is still within bounds
    // with density dependence this should always be the case
    assert(Ni + Ns <= N);
} //end infection_susceptible()

// this will need to be changed if we allow more than 1 plasmid per cell 
void conjugation_infected(int const Receive_idx)
{
    assert(Receive_idx >= 0);
    assert(Receive_idx < Ni);

    // check indeed that both infected
    assert(Infected[Receive_idx].nplasmids > 0);
    assert(Infected[Receive_idx].nplasmids <= nplasmid_max);

    int Donor_idx;
    // if infected has no preference
    // pick a random infected individual 
    // as Donor
    bool haspref = Infected[Receive_idx].p_plasmid || Infected[Receive_idx].p_chr;
    
    if(!haspref)
    {
        std::uniform_int_distribution<int> infected_sampler(0, Ni - 1);
	do Donor_idx = infected_sampler(rng_r);
	while (Donor_idx == Receive_idx);
    }
	    
    // if infected has preference
    // make weighted distribution for 
    // homozygote preference
    // or
    // heterozygote preference
    // and sample from that distribution
    if(haspref)
    {
	if(Infected[Receive_idx].p_plasmid == Infected[Receive_idx].p_chr)	
	{
            std::discrete_distribution <int> infected_attract_dist(attract_homozygote.begin(), attract_homozygote.end());
	    Donor_idx = infected_attract_dist(rng_r);
	}
	else 
	{
            std::discrete_distribution <int> infected_attract_dist(attract_heterozygote.begin(), attract_heterozygote.end());
	    Donor_idx = infected_attract_dist(rng_r);
	}	
    }

    assert(Infected[Donor_idx].nplasmids > 0);
    assert(Infected[Donor_idx].nplasmids <= nplasmid_max);

    // no change in plasmid numbers because 
    // receiver keeps only one plasmid (random choice
    // and donor is giving a copy of plasmid
    // so keeps a copy for themselves 

/* Is there recombination between plasmids before one plasmid chucked out?
*/
    // receiver keeps one plasmid - random choice
    if (uniform(rng_r) < 0.5) 
    	{
    	Infected[Receive_idx].p_plasmid = Infected[Donor_idx].p_plasmid;
    	Infected[Receive_idx].t_plasmid = Infected[Donor_idx].t_plasmid;
	}
    //update attractiveness vectors 
    update_conj_attract(Receive_idx);
    // check whether population is still within bounds
    // with density dependence this should always be the case
    assert(Ni + Ns <= N);
} // end of conjugation_infected()

void loss_plasmid(int const I_idx)
{
    assert(Infected[I_idx].nplasmids > 0);
    assert(I_idx >= 0);
    assert(I_idx < Ni);

    // if we want to have more than 1 plasmid, this bit of code needs the for loop
//    for (int plasmid_idx = 0; plasmid_idx < Infected[I_idx].nplasmids; ++plasmid_idx)
//    {
		    // copy the last plasmid on the stack over the plasmid that's getting lost
		    // then remove the last plasmid on the stack
    Infected[I_idx].plasmid[0] = false;

    --Infected[I_idx].nplasmids;

//    }

	// if no plasmids left, add the individual
	// to the end of the Susceptible stack
	// and write over its place in the infected stack
	// with the last individual in the Infected stack
    if (Infected[I_idx].nplasmids == 0)
    {
        Susceptible[Ns++] = Infected[I_idx];

	//an infected individual is lost
	//update attraction vectors
	update_death_attract(I_idx);

        Infected[I_idx] = Infected[Ni - 1];
    	--Ni;
    }

}// end loss_plasmid()

// recombination between plasmid and chromosome alleles
void recombination(Individual &ind)
{
    assert(ind.nplasmids <= nplasmid_max);
    Individual old_ind = ind;

    if(uniform(rng_r) < r)
    {
        ind.p_plasmid = old_ind.p_chr;
        ind.p_chr = old_ind.p_plasmid;
    }

    if(uniform(rng_r) < r)
    {
        ind.t_plasmid = old_ind.t_chr;
        ind.t_chr = old_ind.t_plasmid;
    }
}// end void recombination()

// birth event of an individual
// is going to be different with 
void birth(Individual &parent, bool parent_susceptible)
{
    Individual kid;
    int idx_kid;

    assert(parent.p_chr == 0 || parent.p_chr ==1);
    assert(parent.t_chr == 0 || parent.t_chr == 1);
    assert(parent.p_plasmid == 0 || parent.p_plasmid ==1);
    assert(parent.t_plasmid == 0 || parent.t_plasmid == 1);

    // mutate the preference and trait alleles
    // in the chromosome
    kid.p_chr = mutation(mu_p, parent.p_chr);
    kid.t_chr = mutation(mu_t, parent.t_chr);

    // for now assuming strict inheritance
    // later on we might want to introduce segregational effects
    kid.nplasmids = 0;

    if (parent_susceptible)
        {
        Susceptible[Ns++] = kid;
        assert(Ns + Ni <= N);
        return; 
        }

    // birth of infected individual
    if (!parent_susceptible) 
	{
	assert(parent.nplasmids > 0);
	for (int plasmid_idx = 0; plasmid_idx < parent.nplasmids; ++plasmid_idx)
	    {
	    // replicate plasmid to offspring
            kid.plasmid[plasmid_idx] = parent.plasmid[plasmid_idx];
	    kid.nplasmids += kid.plasmid[plasmid_idx]; 
	    kid.p_plasmid = parent.p_plasmid;
	    kid.t_plasmid = parent.t_plasmid;
	    recombination(kid);
	    }

        Infected[Ni++] = kid; 
	idx_kid = Ni - 1;
	update_birth_attract(idx_kid); 
        
	//bounds checking
        assert(Ns + Ni <= N);
        assert(kid.p_chr == 0 || kid.p_chr ==1);
        assert(kid.t_chr == 0 || kid.t_chr == 1);
        assert(kid.p_plasmid == 0 || kid.p_plasmid ==1);
        assert(kid.t_plasmid == 0 || kid.t_plasmid == 1);
        return;
	} // end if (!parent_susceptible)

} // end of birth()

// death of susceptible individual
void death_susceptible(int const random_susceptible)
{

    assert(random_susceptible >= 0);
    assert(random_susceptible < Ns);
    Susceptible[random_susceptible] = Susceptible[Ns - 1];
    --Ns;

    assert(Ns >= 0);
    assert(Ns + Ni <= N);
}// end death_susceptible()

// in this model the death rate for susceptible 
// and infected individuals is the same
void death_infected(int const I_idx)
{
    assert(Infected[I_idx].nplasmids > 0);
    assert(I_idx >= 0);
    assert(I_idx < Ni);
    update_death_attract(I_idx);
    Infected[I_idx] = Infected[Ni - 1];
    --Ni;

    assert(Ni >= 0);
    assert(Ns + Ni <= N);
}

/* Leaving this here for now
 * we will need it when we allow 
 * for more than one plasmid to co-exist in a cell
 */


// death of an infected individual at location I_idx

// write headers to the datafile
void write_data_headers(std::ofstream &data_file)
{
    data_file << "time;Ns;Ni;freq_p2_all;freq_p2_infected;freq_p2_plasmid;freq_t2_all;freq_t2_infected;freq_t2_plasmid;var_p2_all;var_p2_infected;var_p2_plasmid;var_t2_all;var_t2_infected;var_t2_plasmid;mean_nplasmid;var_nplasmid;" << std::endl;
} // end write_data_headers()

// fecundity function that accounts for costly trait and preference,
// as in Gandon & Vale eq (2)
double b_Susceptible(bool const trait, bool const pref)
{
    return(bmax * exp(- c * pref - epsilon * trait));
}

double b_Infected(bool const trait_pl, bool const pref_pl, bool const trait_chr, bool const pref_chr)
{
   // check if individual has trait
   // and preference
    bool hastrait = trait_pl || trait_chr;
    bool haspref = pref_pl || pref_chr;
   
   // should the trait or preference be 
   // multiplied by global variables for dominance coefficient? 
    double l_calc; // trait dominance coeff 
    double h_calc; // pref dominance coeff 

    if(trait_pl != trait_chr)
        {
	    l_calc = l;
        }
     else l_calc = 1;

    if(pref_pl != pref_chr)
        {
	    h_calc = h;
        }
    else h_calc = 1;

    return(bmax * exp(- c * h_calc * haspref - epsilon * l_calc * hastrait - delta));
} // end of b_Infected()

// calculate the odds of a plasmid loss event
// and change loss_rate accordingly
// it will be equal to gamma while there's only 1 plasmid / individual
void gamma_loss(double &loss_rate, Individual &host)
{
    assert(host.nplasmids <= nplasmid_max);
    assert(host.nplasmids > 0);

    loss_rate = host.nplasmids * loss_gamma;
}

// setup a distribution of events and choose
// which event to do
void event_chooser(int const time_step)
{
    // some auxiliary variables to store temporary values about rates
    double rate_birth, rate_loss, rate_infect, conjugation_rate, death_rate_S, death_rate_I ;

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

    // declare vectors containing rates of
    // plasmid loss 
    std::vector <double> loss_rates;

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
    // go through all susceptibles and calculate birth
    // and infection rates
    for (int S_idx = 0; S_idx < Ns; ++S_idx)
    {
/* 
        if (Susceptible[S_idx].x < 0.0)
        {
            std::cout << "time: " << time_step << " Ns: " << Ns << " inf_idx: " << S_idx << " " << Susceptible[S_idx].x << std::endl;
        }
*/

        assert(Susceptible[S_idx].nplasmids == 0);

        // 0. Birth of susceptibles
        rate_birth = b_Susceptible(Susceptible[S_idx].t_chr, Susceptible[S_idx].p_chr) * (1.0 - kappa * (Ns + Ni));
        birth_rates_susceptible.push_back(rate_birth);

        total_rates[0] += rate_birth;

        
    } // end for S_idx

    // 1. infection of susceptible by plasmid
    // in the fisherian sexsel model 'infection' occurs via conjugation
    // not environmental uptake
    // hence the rate of infections should depend on number of 
    // bacteria already infected
    // COME BACK to this spot after reading Keeling & Rohani
    rate_infect =  ((Ns*Ni)/ N) * psi;
        
    total_rates[1] = rate_infect;

    // now go through the infected hosts 
    for (int inf_idx = 0; inf_idx < Ni; ++inf_idx)
    {
        // bounds checking. If there is a buffer overflow
        // these numbers typically are not between bounds
        //std::cout << "time: " << time_step << " Ni: " << Ni << " inf_idx: " << inf_idx << " nplasmids: " << Infected[inf_idx].nplasmids << std::endl;
        assert(Infected[inf_idx].nplasmids > 0.0);
        assert(Infected[inf_idx].nplasmids <= nplasmid_max);
        
   /*     if (Infected[inf_idx].x < 0.0)
        {
            std::cout << "time: " << time_step << " Ni: " << Ni << " inf_idx: " << inf_idx << " " << Infected[inf_idx].x << std::endl;
        }
*/

        // 2. birth infected host   
        rate_birth = b_Infected(Infected[inf_idx].t_plasmid, Infected[inf_idx].p_plasmid, Infected[inf_idx].t_chr, Infected[inf_idx].p_chr) * (1.0 - kappa * (Ns + Ni));
        
        birth_rates_infected.push_back(rate_birth);

        total_rates[2] += rate_birth;

        // 3.  loss of a plasmid 
        // a single plasmid loss event
        gamma_loss(rate_loss, Infected[inf_idx]);

	//update individual loss rates
	// as well as cumulative loss rate
        loss_rates.push_back(rate_loss);
        total_rates[3] += rate_loss;

    }

    // 4. conjugation between infected
    conjugation_rate = ((Ni*Ni)/N) * psi;
        
    total_rates[4] = conjugation_rate;

    // 5. Deaths susceptibles
    death_rate_S = Ns * d;

    total_rates[5] = death_rate_S;

    // 6. death rate infected
    death_rate_I = Ni * d; 

    total_rates[6] = death_rate_I;

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
                std::cout << "Birth of susceptible" << std::endl;
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

		// note that S_idx is an index from 0 to Ns - 1
		// thus it does not include Ns itself
		assert(S_idx >= 0);
		assert(S_idx < Ns);

		// execute the birth() function which also updates the stats
		birth(Susceptible[S_idx], true);

		break;
            }

        case 1: // infection of a susceptible 
            {
                std::cout << "Infection of susceptible" << std::endl;
                // choose a susceptible individual at random 
                // to get infected
                // then draw from the probability distribution
                // to determine the individual that gets infected
                std::uniform_int_distribution<int> susceptible_sampler(0, Ns - 1);
                S_idx = susceptible_sampler(rng_r);

                // bounds checking
                assert(S_idx >= 0);
                assert(S_idx < Ns);

                infection_susceptible(S_idx);

                break;
            }

        case 2: // birth infected host
            {
                std::cout << "Birth of infected" << std::endl;
		assert(birth_rates_infected.size() == Ni);
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

        case 3: // loss of a single plasmid
            {
                std::cout << "Loss of plasmid" << std::endl;
                // set up probability distribution that determines
                // which individual will lose a plasmid
                std::discrete_distribution <int> loss_plasmid_dist(
                        loss_rates.begin()
                        ,loss_rates.end());

                // draw an individual from that distribution which is going
                // to lose a plasmid
                I_idx = loss_plasmid_dist(rng_r);

                // bounds checking
                assert(I_idx >= 0);
                assert(I_idx < Ni);
               
                // perform the actual plasmid loss 
                loss_plasmid(I_idx);
            
                break;
            }

        case 4: // conjugation infected 
            {
                std::cout << "Conjugation of infected" << std::endl;
                //randomly pick
                //which individual will be a receiver 
                // draw the individual that is going to get co-infected
                std::uniform_int_distribution<int> infected_sampler(0, Ni - 1);
                I_idx = infected_sampler(rng_r);

                // bounds check
                assert(I_idx >= 0);
                assert(I_idx < Ni);

                // perform the actual conjugation
		// where a donor will be selected 
		// according to the receiver's preference
                conjugation_infected(I_idx);

                break;
            }

        case 5: // death susceptible
	    {
            // each susceptible has the same chance of dying
            // so we do not need to set up a probability distribution
            // determining which susceptible is more likely to die relative
            // to others, we simply pick a 
            // random individual
                std::cout << "Death susceptible" << std::endl;

                std::uniform_int_distribution<int> susceptible_sampler(0, Ns - 1);
                S_idx = susceptible_sampler(rng_r);
                assert(S_idx >= 0);
                assert(S_idx < Ns);
                
                death_susceptible(S_idx);

                break;
	    }

        case 6:// death infected
            {
                std::cout << "Death infected" << std::endl;
                //currently, infected individuals
		//all have same chance of dying
		//so pick one at random 
                std::uniform_int_distribution<int> infected_sampler(0, Ni - 1);
                I_idx = infected_sampler(rng_r);

                // bounds check
                assert(I_idx >= 0);
                assert(I_idx < Ni);
                
                // perform actual death 
                death_infected(I_idx);

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

    double mean_n_plasmid = 0.0;
    double ss_n_plasmid = 0.0;

    double mean_freq_p2_total = 0.0;
    double ss_freq_p2_total = 0.0;

    double mean_freq_p2_infected = 0.0;
    double ss_freq_p2_infected = 0.0;

    double mean_freq_p2_plasmid = 0.0;
    double ss_freq_p2_plasmid = 0.0;

    double mean_freq_t2_total = 0.0;
    double ss_freq_t2_total = 0.0;

    double mean_freq_t2_infected = 0.0;
    double ss_freq_t2_infected = 0.0;

    double mean_freq_t2_plasmid = 0.0;
    double ss_freq_t2_plasmid = 0.0;

    double p2, p2_plasmid, t2, t2_plasmid, n_plasmid;

    for (int I_idx = 0; I_idx < Ni; ++I_idx) {
	
	assert(Infected[I_idx].p_plasmid <= 1);
	assert(Infected[I_idx].p_chr <= 1);
	assert(Infected[I_idx].t_plasmid <= 1);
	assert(Infected[I_idx].t_chr <= 1);

        p2 = Infected[I_idx].p_plasmid + Infected[I_idx].p_chr;
        p2_plasmid = Infected[I_idx].p_plasmid;

	assert(p2 <= 2);
	assert(p2_plasmid <=1);

//	std::cout << "I_idx " << I_idx << "allele p_plasmid  " <<  Infected[I_idx].p_plasmid << std::endl;
//	std::cout << "I_idx " << I_idx << "p2_plasmid  " <<  p2_plasmid << std::endl;


        mean_freq_p2_total += p2;
        ss_freq_p2_total += p2 * p2;

        mean_freq_p2_infected += p2;
        ss_freq_p2_infected += p2 * p2;

        mean_freq_p2_plasmid += p2_plasmid;
//	std::cout << "mean frequency p2 plasmid sum " << mean_freq_p2_plasmid << std::endl;
        ss_freq_p2_plasmid += p2_plasmid * p2_plasmid;

        t2 = Infected[I_idx].t_plasmid + Infected[I_idx].t_chr;
        t2_plasmid = Infected[I_idx].t_plasmid;

	assert(t2 <=2);
	assert(t2_plasmid <=1);

        mean_freq_t2_total += t2;
        ss_freq_t2_total += t2 * t2;

        mean_freq_t2_infected += t2;
        ss_freq_t2_infected += t2 * t2;

        mean_freq_t2_plasmid += t2_plasmid;
        ss_freq_t2_plasmid += t2_plasmid * t2_plasmid;

	n_plasmid = Infected[I_idx].nplasmids;
        mean_n_plasmid += n_plasmid;
        ss_n_plasmid += n_plasmid * n_plasmid;
    }

    assert(mean_freq_p2_infected <= 2 * Ni);
    assert(mean_freq_t2_infected <= 2 * Ni);
    assert(mean_freq_p2_plasmid <=  Ni);
    assert(mean_freq_t2_plasmid <=  Ni);
    assert(mean_n_plasmid <= Ni);
    
    mean_freq_p2_infected /= 2*Ni;
    mean_freq_p2_plasmid /= Ni;

    mean_freq_t2_infected /= 2*Ni;
    mean_freq_t2_plasmid /= Ni;
    
    mean_n_plasmid /= Ni;

    double var_p2_infected = ss_freq_p2_infected / 2*Ni 
        - mean_freq_p2_infected * mean_freq_p2_infected;

    double var_t2_infected = ss_freq_t2_infected / 2*Ni 
        - mean_freq_t2_infected * mean_freq_t2_infected;

    double var_p2_plasmid = ss_freq_p2_plasmid / Ni 
        - mean_freq_p2_plasmid * mean_freq_p2_plasmid;

    double var_t2_plasmid = ss_freq_t2_plasmid / Ni 
        - mean_freq_t2_plasmid * mean_freq_t2_plasmid;

    double var_n_plasmid = ss_n_plasmid / Ni 
	    - mean_n_plasmid * mean_n_plasmid;

    for (int S_idx = 0; S_idx < Ns; ++S_idx) {

	assert(Infected[S_idx].p_chr <= 1);

	assert(Infected[S_idx].t_chr <= 1);

        p2 = Susceptible[S_idx].p_chr;

        mean_freq_p2_total += p2;
        ss_freq_p2_total += p2 * p2;

        t2 = Susceptible[S_idx].t_chr;

        mean_freq_t2_total += t2;
        ss_freq_t2_total += t2 * t2;

    }
    
    assert(mean_freq_p2_total <= 2 * Ni + Ns); 
    assert(mean_freq_t2_total <= 2 * Ni + Ns); 

    mean_freq_p2_total /= (2*Ni + Ns);
    mean_freq_t2_total /= (2*Ni + Ns);

    double var_p2_total = ss_freq_p2_total / (2*Ni + Ns) 
        - mean_freq_p2_total * mean_freq_p2_total;

    double var_t2_total = ss_freq_t2_total / (2*Ni + Ns)
    	    - mean_freq_t2_total * mean_freq_t2_total;

    data_file << time_step << ";"
        << Ns << ";"
        << Ni << ";"
	<< mean_freq_p2_total << ";"
	<< mean_freq_p2_infected << ";"
	<< mean_freq_p2_plasmid << ";"
	<< mean_freq_t2_total << ";"
	<< mean_freq_t2_infected << ";"
	<< mean_freq_t2_plasmid << ";"
	<< var_p2_total << ";"
	<< var_p2_infected << ";"
	<< var_p2_plasmid << ";"
	<< var_t2_total << ";"
	<< var_t2_infected << ";"
	<< var_t2_plasmid << ";"
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



