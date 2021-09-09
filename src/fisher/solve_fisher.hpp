#ifndef ITERATE_FISHER_HPP
#define ITERATE_FISHER_HPP

#include <string>
#include <iostream>
#include <fstream>


enum genotype {
    t1p1 = 0,
    t2p1 = 1,
    t1p2 = 2,
    t2p2 = 3
};


// evolution of the following chromosomal
// genotypes:
// (p1, t1), (p1, t2), (p2, t2)
class SolveFisher
{
    private:
        // vanish boundary below which we consider that there
        // is no change in x(t+1) - x(t)
        static constexpr double vanish_bound{1e-07};
        static constexpr double eul{0.01};

        void init_arguments(int argc, char ** argv);

        // allocate vectors S for the 
        // three chromosomal genotypes
        double S[4];

        // allocate vector I for the 
        // 3 chromosomal x 3 plasmid genotypes
        double I[4][4];

        // fecundities for each of the 4 genos
        double bS[4];
        double bI[4][4];

        // infection rates with S recipient and I donor
        double beta_SxI[4][4][4];

        // infection rates with I recipient and I donor
        double beta_IxI[4][4][4][4];

        long int max_time;
        double kappa;
        double gamma;
        double d;
        double delta;
        double ht;
        double hp;
        double r;
        double a;
        double mu_t[2];
        double mu_p[2];

        // infection prob
        double pi;

        double N;

        // cost orn, pref
        double cp;
        double ct;

        // initial freq of p2 and t2
        double p2_t0;
        double t2_t0;

        // initial freq of linkage disequilibrium
        double D_t0;

        int skip_output;

        std::string base_name;

        // array to translate boolean allelic values to genotype
        // first index: ti
        // second index: pi
        genotype allele2genotypes[2][2];

        bool has_p2[4];
        bool has_t2[4];

        void write_parameters(std::ofstream &data_file);
        void write_data(std::ofstream &data_file);

    public :
        SolveFisher(int argc, char **argv); // c'tor

        // solve the system of differential equations
        void solveSys();
};



#endif
