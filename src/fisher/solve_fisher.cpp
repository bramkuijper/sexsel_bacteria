#include <iostream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include "solve_fisher.hpp"


void SolveFisher::SolveFisher(int argc, char **argv) :
    S{0.0,0.0,0.0}
    ,I{{0.0,0.0,0.0,0.0}
        ,{0.0,0.0,0.0,0.0}
        ,{0.0,0.0,0.0,0.0}
        ,{0.0,0.0,0.0,0.0}}
    ,bS{0.0,0.0,0.0}
    ,bI{{0.0,0.0,0.0,0.0}
        ,{0.0,0.0,0.0,0.0}
        ,{0.0,0.0,0.0,0.0}
        ,{0.0,0.0,0.0,0.0}}
    ,beta_SxI{}
    ,beta_IxI{}
    ,kappa{0.0}
    ,delta{0.0}
    ,max_time{0.0}
    ,gamma{0.0}
    ,d{0.0}
    ,pi{0.0}
    ,cp{0.0}
    ,ct{0.0}
    ,ht{0.0}
    ,hp{0.0}
    ,mu_t{0.0,0.0}
    ,mu_p{0.0,0.0}
    ,a{0.0}
    ,N{0.0}
    ,p2_t0{0.0} 
    ,t2_t0{0.0} 
    ,D_t0{0.0} 
    ,base_name{} 
    ,allele2genotypes{{t1p1,t1p2},{t2p1,t2p2}}
    ,has_p2{false,false,true,true}
    ,has_t2{false,true,false,true}
{
    init_arguments(argc, argv);
}

// initialize arguments from the command line
// and fill vectors that are functions of basic 
// parameters
void SolveFisher::init_arguments(int argc, char **argv)
{
    kappa = atof(argv[1]);
    gamma = atof(argv[1]);
    delta = atof(argv[1]);
    max_time = atof(argv[1]);
    d = atof(argv[1]);
    N = atof(argv[2]);
    p2_t0 = atof(argv[3]);
    t2_t0 = atof(argv[3]);
    D_t0 = atof(argv[3]);
    cp = atof(argv[3]);
    ct = atof(argv[3]);
    pi = atof(argv[3]);
    a = atof(argv[3]);
    hp = atof(argv[3]);
    ht = atof(argv[3]);
    mu_t[0] = atof(argv[3]);
    mu_t[1] = atof(argv[3]);
    mu_p[0] = atof(argv[3]);
    mu_p[1] = atof(argv[3]);

    base_name = argv[10];

    // fill vectors for parameters that will not change
    // during the iteration, i.e., fecundity and transmission rates
    for (int geno_recip_chr_idx = 0; 
            geno_recip_chr_idx < 4; 
                ++geno_recip_chr_idx)
    {
        // first calculate fecundity of susceptibles

        // get the genotype number
        genotype geno = static_cast<genotype>(geno_recip_chr_idx);

        // number of t2s on recipient's chromosome (0 or 1)
        int nt2_recip_chr = geno == p1t2 || geno == t2p2;

        // number of t2s on recipient's chromosome (0 or 1)
        int np2_recip_chr = geno == t1p2 || geno == t2p2;

        // fecundity of susceptible depends on number of t2s and p2s
        bS[geno_recip_chr_idx] = exp(-cp*np2_recip_chr -ct*nt2_recip_chr);
        
        // second, calculate fecundity of infected individuals

        // now go through the different plasmid genotypes of the donor
        // so that we can calculate fecundities of infecteds
        for (int geno_donor_plm_idx = 0; geno_donor_plm_idx < 4; 
                ++geno_donor_plm_idx)
        {
            // get plasmid genotype of donor
            genotype geno_plasmid_donor = 
                static_cast<genotype>(geno_donor_plm_idx);

            // count number of t2s on plasmid of donor
            int nt2_donor_plm = geno_plasmid_donor == p1t2 
                || geno_plasmid_donor == t2p2;

            // count number of p2s on plasmid of donor
            int np2_donor_plm = geno_plasmid_donor == t1p2 
                || geno_plasmid_donor == t2p2;

            // fecundity of infected
            bI[geno_recip_chr_idx][geno_donor_plm_idx] = 
                exp(-cp*(np2_recip_chr+np2p) 
                        -ct*(nt2_recip_chr+nt2_donor_plm));

            // now loop through chromosomal genotypes of infected
            // and calculate infection rates
            for (int geno_donor_chr_idx = 0; 
                    geno_donor_chr_idx < 4; ++geno_donor_chr_idx)
            {
                genotype geno_donor_chr = 
                    static_cast<genotype>(geno_donor_chr_idx);

                int nt2_donor_chr = geno_donor_chr == p1t2 || 
                    geno_donor_chr == t2p2;

                int nt2_donor_chr = geno_donor_chr == t1p2 || 
                    geno_donor_chr == t2p2;

                // by default set the infection rate to 0
                beta_SxI[geno_recip_chr_idx][
                    geno_donor_plm_idx][geno_donor_chr_idx] = 1.0;

                // if recipient is choosy on either chrom or plasmid
                // and something in donor that looks like an orament
                if (np2_recip_chr > 0 && nt2_donor_plm + nt2_donor_chr > 0)
                {
                    beta_SxI[geno_recip_chr_idx][geno_donor_plm_idx][geno_donor_chr_idx] += 
                        nt2_donor_plm + (n_t2_donor chr == 1.0 ? 
                            ht * a // donor heterozygous t1 t2
                            :
                            a); // donor homozygous
                }

                for (int geno_recip_plm_idx = 0; 
                        geno_recip_plm_idx < 4; ++geno_recip_plm_idx)
                {
                    genotype geno_recip_plm = 
                        static_cast<genotype>(geno_recip_plm_idx);

                    int nt2_recip_plm = geno_recip_plm == t2p1 || geno_recip_plm == t2p2;
                    int np2_recip_plm = geno_recip_plm == t1p2 || geno_recip_plm == t2p2;

                    // infection rate from the perspective of an infected recipient
                    // getting infected by an infected donor
                    beta_IxI[geno_recip_chr_idx][geno_recip_plm_idx][
                        geno_donor_chr_idx][geno_donor_plm_idx] = 1.0;

                    // recipient heterozygous for p2
                    // any ornamented donor
                    if (np2_recip_chr + np2_recip_plm == 1 
                            && nt2_donor_plm + nt2_donor_chr > 0)
                    {
                        beta_IxI[geno_recip_chr_idx][geno_recip_plm_idx][
                            geno_donor_chr_idx][geno_donor_plm_idx] +=
                               nt2_donor_plm + nt2_donor_chr == 1 ?
                                ht * hp * a
                                :
                                hp * a;
                    }
                    else if (np2_recip_chr + np2_recip_plm == 2 && 
                            nt2_donor_plm + nt2_donor_chr > 0)
                    {
                        beta_IxI[geno_recip_chr_idx][geno_recip_plm_idx][
                            geno_donor_chr_idx][geno_donor_plm_idx] +=
                               nt2_donor_plm + nt2_donor_chr == 1 ?
                                ht * a
                                :
                                a;
                    }
                } // end for int geno_recip_plm_idx
            } // end for geno_donor_chr_idx
        } // end for geno_donor_plm_idx
    } // end for geno_recip_chr_idx

    assert(beta_IxI[3][3][3][3] == a);
    assert(beta_IxI[3][0][3][3] == 1.0);

}//end init_arguments

// write all the parameters to the file data_file
void SolveFisher::write_parameters(std::ofstream &data_file)
{
    for (int geno_chr_idx = 0; geno_chr_idx < 4; ++geno_chr_idx)
    {
        genotype geno_chr = static_cast<genotype>(geno_chr_idx);
        std::string tval_chr = geno_chr == t2p1 || geno_chr == t2p2 ? "t2" : "t1";
        std::string pval_chr = geno_chr == t1p2 || geno_chr == t2p2 ? "p2" : "p1";

        data_file << "fec_S" 
            << tval_chr 
            << pval_chr 
            << ";" 
            << bS[geno_chr_idx] << std::endl;

        for (int geno_plm_idx = 0; geno_plm_idx < 4; ++geno_plm_idx)
        {
            genotype geno_plm = static_cast<genotype>(geno_plm_idx);
            std::string tval_plm = geno_plm == t2p1 || geno_plm == t2p2 ? "t2" : "t1";
            std::string pval_plm = geno_plm == t1p2 || geno_plm == t2p2 ? "p2" : "p1";

            data_file << "fec_I_" 
                << tval_chr 
                << pval_chr 
                << tval_plm 
                << pval_plm << ";" << 
                bI[geno_chr][geno_plm_idx] << std::endl;
        
            for (int geno_chr_idx2 = 0; 
                    geno_chr_idx2 < 4; ++geno_chr_idx2)
            {
                genotype geno_donor_chr = static_cast<genotype>(geno_chr_idx2);

                std::string tval_donor_chr = geno_donor_chr == t2p1 || geno_donor_chr == t2p2 ? 
                    "t2" : "t1";
                std::string pval_donor_chr = geno_donor_chr == t1p2 || geno_donor_chr == t2p2 ? 
                    "p2" : "p1";

                data_file << "beta_SxI_" 
                    << tval_chr 
                    << pval_chr << "_"
                    << tval_donor_chr 
                    << pval_donor_chr 
                    << tval_plm 
                    << pval_plm << ";" 
                    << beta_SxI[geno_chr_idx][geno_chr_idx2][geno_plm_idx]
                    << std::endl;

                for (int geno_plm_idx2 = 0; 
                        geno_plm_idx2 < 4; ++geno_plm_idx2)
                {
                    genotype geno_recip_plm = static_cast<genotype>(geno_plm_idx2);

                    std::string tval_recip_plm = 
                        geno_recip_plm == t2p1 || geno_recip_plm == t2p2 ? 
                        "t2" : "t1";

                    std::string pval_recip_plm = 
                        geno_recip_plm == t1p2 || geno_recip_plm == t2p2 ? 
                        "p2" : "p1";

                    data_file << "beta_IxI_" 
                        << tval_chr 
                        << pval_chr 
                        << tval_recip_plm
                        << pval_recip_plm
                        << "_"
                        << tval_donor_chr 
                        << pval_donor_chr 
                        << tval_plm 
                        << pval_plm << ";" 
                        << beta_IxI[geno_chr_idx][geno_plm_idx2][geno_chr_idx2][geno_plm_idx]
                        << std::endl;

                } // end for (int geno_plm_idx2
            } // end for (int geno_chr_idx2 = 0
        } // end for (int geno_plm_idx = 0
    } // end for (int geno_chr_idx

    data_file << "kappa;" << kappa << std::endl
        << "gamma;" << gamma << std::endl
        << "d;" << d << std::endl
        << "ht;" << ht << std::endl
        << "hp;" << hp << std::endl
        << "a;" << a << std::endl
        << "mu_t1;" << mu_t[0] << std::endl
        << "mu_t2;" << mu_t[1] << std::endl
        << "mu_p1;" << mu_p[0] << std::endl
        << "mu_p2;" << mu_p[1] << std::endl
        << "pi;" << pi << std::endl
        << "N;" << N << std::endl
        << "cp;" << cp << std::endl
        << "ct;" << ct << std::endl
        << "p2_t0;" << p2_t0 << std::endl
        << "t2_t0;" << t2_t0 << std::endl
        << "D_t0;" << D_t0 << std::endl;
}// end write_parameters

// attempt to numerically solve the 
// system of differential equations
void SolveFisher::solveSys()
{
    double dSdt[4];
    double dIdt[4][4];

    bool converged;

    std::ofstream output_file(base_name);

    // write parameters to output file
    write_parameters(output_file);


    for (int time_idx = 0; time_idx < max_time; ++time_idx)
    {
        // calculate rates of plasmid loss
        double sum_plasmid_loss[4] = {0.0,0.0,0.0,0.0};
        
        // calculate total force of infection
        // either in terms of loss of a susceptible of genotype i
        double total_force_of_infection_lossS[4] = {0.0,0.0,0.0,0.0};

        // ... or rate of gain of new infecteds of 
        // genotype i and plasmid genotype j
        double total_force_of_infection_gainI[4][4];

        // total force of infection between two infected individuals
        double total_force_of_infection_IxI[4][4][4][4];

        // recombination between chromosome and plasmid in infected individuals;
        double recombination[4][4];

        // reset population size (N) 
        // count and update it
        N = 0;

        // - update value for N
        // - calculate plasmid loss rate
        // - reset force of infection calculation
        for (int genotype_chr_idx = 0; genotype_chr_idx < 4; ++genotype_chr_idx)
        {
            N += S[genotype_chr_idx];
            
            total_force_of_infection_lossS[genotype_chr_idx] = 0.0;

            // level 2: donor plasmid 
            for (int genotype_plm_idx = 0; 
                    genotype_plm_idx < 4; ++genotype_plm_idx)
            {
                // update N with infecteds
                N+=I[genotype_chr_idx][genotype_plm_idx];

                // set total rate of gain in I due to infections to 0 and
                // calculate value in the next loop 
                total_force_of_infection_gainI[genotype_chr_idx][genotype_plm_idx] = 0.0;

                sum_plasmid_loss[genotype_chr_idx] += 
                    gamma * I[genotype_chr_idx][genotype_plm_idx];
            }
            assert(sum_plasmid_loss[genotype_chr_idx] >= 0);
        }

        assert(N > 0);

        // level 1 recipient chromosome
        for (int genotype_chr_idx = 0; 
                genotype_chr_idx < 4; ++genotype_chr_idx)
        {
            // level 2: donor chromosome
            for (int genotypeI_chromosome_idx = 0; 
                    genotypeI_chromosome_idx < 4; ++genotypeI_chromosome_idx)
            {
                // level 3: donor plasmid
                for (int genotypeI_plasmid_idx = 0; 
                        genotypeI_plasmid_idx < 4; ++genotypeI_plasmid_idx)
                {
                    // calculate gain_rate of I due to infection of susceptibles
                    total_force_of_infection_gainI[
                        genotypeI_chromosome_idx][genotypeI_plasmid_idx] += 
                            beta_SxI[genotype_chr_idx][
                            genotypeI_chromosome_idx][genotypeI_plasmid_idx] * 
                                S[genotype_chr_idx] / N;

                    // calculate loss rate of S due to infection of susceptibles
                    total_force_of_infection_lossS[genotype_chr_idx] += 
                        beta_SxI[genotype_chr_idx][
                            genotypeI_chromosome_idx][genotypeI_plasmid_idx
                            ] * I[genotypeI_chromosome_idx][genotypeI_plasmid_idx] / N;
                }

            }

        }

        // then calculate the infected x infected rate
        // for sake of clarity I do this in separate loop

        // level 1 recipient chromosome
        for (int geno_recip_chr_idx = 0; geno_recip_chr_idx < 4; ++geno_recip_chr_idx)
        {
            // level 2 recipient plasmid
            for (int geno_recip_plm_idx = 0; geno_recip_plm_idx < 4; ++geno_recip_plm_idx)
            {
                bool p2_chr = has_p2[geno_recip_chr_idx];
                bool p2_plm = has_p2[geno_recip_plm_idx];
                bool t2_chr = has_t2[geno_recip_chr_idx];
                bool t2_plm = has_t2[geno_recip_plm_idx];

                // only recombination if alleles on both chromosomes are different
                if (p2_chr != p2_plm && t2_chr != t2_plm)
                {
                    // two types of influx
                    // if this is t1p1 x t2p2 we have
                    // t2p1 and t1p2, each with freq 1/2 * r
                    genotype recombinant1 = allele2genotypes[p2_chr][t2_plm];
                    genotype recombinant2 = allele2genotypes[p2_plm][t2_chr];


                    recombination_in[geno_recip_chr_idx][geno_recip_plm_idx] +=
                        0.5 * r * recombination[geno_recip_chr_idx][geno_recip_plm_idx];

                    recombination_out[geno_recip_chr_idx][geno_recip_plm_idx] -=

                }


                // level 3: donor chromosome
                for (int geno_donor_chr_idx = 0; 
                        geno_donor_chr_idx < 4; ++geno_donor_chr_idx)
                {
                    // level 4: donor plasmid
                    for (int geno_donor_plm_idx = 0; 
                            geno_donor_plm_idx < 4; ++geno_donor_plm_idx)
                    {
                        total_force_of_infection_IxI[
                            geno_recip_chr_idx][
                                geno_recip_plm_idx][
                                    geno_donor_chr_idx][
                                        geno_donor_plm_idx] = 
                            beta_IxI[
                                geno_recip_chr_idx][
                                    geno_recip_plm_idx][
                                        geno_donor_chr_idx][
                                            geno_donor_plm_idx] *
                            I[geno_recip_chr_idx][
                                    geno_recip_plm_idx] *
                                I[geno_donor_chr_idx][
                                    geno_donor_plm_idx] / N;
                    }
                }
            }
        }

        double total_mutations = 0.0;

        // go through the 4 genotypes on the autosome:
        // (p1, t1), (p1, t2), (p2,t1), (p2, t2)
        for (int genotype_idx = 0; genotype_idx < 4; ++genotype_idx)
        {
            // get the current genotype as 
            // the corresponding enum value
            genotype geno_chr = static_cast<genotype>(genotype_idx);

            // genotype contains t2?
            bool is_t2 = geno_chr == t2p1 || geno_chr == t2p2;

            // genotype contains p2?
            bool is_p2 = geno_chr == t1p2 || geno_chr == t2p2;

            // genotype same as current but opposite t allele
            // (we need to have this to obtain mutation rates)
            genotype geno_chr_oppT = allele2genotypes[!is_t2][is_p2];
            
            // genotype same as current but opposite p allele
            // (we need to have this to obtain mutation rates)
            genotype geno_chr_oppP = allele2genotypes[is_t2][!is_p2];

            // check that total sum of mutations is indeed zero
            // (for the sake of debugging)
            total_mutations += mu_t[!is_t2] * S[geno_chr_oppP] + mu_p[!is_p2] * S[geno_chr_oppP];
            total_mutations -= (mu_t[is_t2] + mu_p[is_p2]) * S[genotype_idx];

            // dSdt
            dSdt[genotype_idx] = 
                // 1. birth of new susceptibles
                bS[genotype_idx] * (1.0 - kappa * N) * S[genotype_idx] +

                // 2. gain in susceptibles due to loss of plasmid
                sum_plasmid_loss[genotype_idx] +

                // 3. mutation incoming
                + mu_t[!is_t2] * S[geno_chr_oppP] + mu_p[!is_p2] * S[geno_chr_oppP]

                // 4. mutation outgoing
                - (mu_t[is_t2] + mu_p[is_p2]) * S[genotype_idx] 
                
                // 5. loss due to deaths
                - d * S[genotype_idx] 
                
                // 6. loss due to infections
                - (1-pi) * total_force_of_infection_lossS[genotype_idx] * S[genotype_idx]; 

            for (int plasmid_idx = 0; plasmid_idx < 3; ++plasmid_idx)
            {
                genotype geno_plm = static_cast<genotype>(plasmid_idx);

                bool is_t2_plm = geno_plm == t2p1 || geno_plm == t2p2;
                bool is_p2_plm = geno_plm == t1p2 || geno_plm == t2p2;

                // genotype same as current but for t allele
                genotype geno_plm_oppT = allele2genotypes[!is_t2_plm][is_p2_plm];
                // genotype same as current but for p allele
                genotype geno_plm_oppP = allele2genotypes[is_t2_plm][!is_p2_plm];

                total_mutations +=

                    // 2a. influx due mutation of t on chromosome
                    mu_t[!is_t2] * I[geno_chr_oppT][plasmid_idx]

                    // 2b. influx due mutation of p on chromosome
                    + mu_p[!is_p2] * I[geno_chr_oppP][plasmid_idx]

                    // 2c. influx due to mutation of t on plasmid
                    + mu_t[!is_t2_plm] * I[genotype_idx][geno_plm_oppT]

                    // 2d. influx due to mutation of p on plasmid
                    + mu_p[!is_p2_plm] * I[genotype_idx][geno_plm_oppP]

                    // 2e. outflux due to mutation on t,p on chrom or plasmid
                    - (mu_t[is_t2] + mu_t[is_t2_plm] + mu_p[is_p2] + mu_p[is_p2_plm]) *
                        I[genotype_idx][plasmid_idx];


                dIdt[genotype_idx][plasmid_idx] = 
                    // 1. birth
                    bI[genotype_idx][plasmid_idx] * (1.0 - kappa * N) * 
                        I[genotype_idx][plasmid_idx] 

                    // 2a. influx due mutation of t on chromosome
                    + mu_t[!is_t2] * I[geno_chr_oppT][plasmid_idx]

                    // 2b. influx due mutation of p on chromosome
                    + mu_p[!is_p2] * I[geno_chr_oppP][plasmid_idx]

                    // 2c. influx due to mutation of t on plasmid
                    + mu_t[!is_t2_plm] * I[genotype_idx][geno_plm_oppT]

                    // 2d. influx due to mutation of p on plasmid
                    + mu_p[!is_p2_plm] * I[genotype_idx][geno_plm_oppP]

                    // 2e. outflux due to mutation on t,p on chrom or plasmid
                    - (mu_t[is_t2] + mu_t[is_t2_plm] + mu_p[is_p2] + mu_p[is_p2_plm]) *
                        I[genotype_idx][plasmid_idx]

                    // 3. gain in infecteds from susceptible
                    + (1-pi) * total_force_of_infection_gainI[genotype_idx][plasmid_idx] *
                        I[genotype_idx][plasmid_idx]

                    // 4. loss in infecteds
                    - (d + gamma) * I[genotype_idx][plasmid_idx];

            } // end for plasmid_idx

        } // end for int genotype_idx

        converged = true;

        for (int genotype_idx = 0; genotype_idx < 4; ++genotype_idx)
        {
            // S[i](t+1) = S[i] + eul * dSdt[i]
            // S[i](t+1) - S[i] < vanish_bound
            // eul*dSdt < vanish_bound
            if (fabs(eul*dSdt[genotype_idx]) > vanish_bound)
            {
                converged = false;
            }

            // update value of S[i]
            S[genotype_idx] += eul*dSdt[genotype_idx];

            for (int plasmid_idx = 0; plasmid_idx < 4; ++plasmid_idx)
            {
                if (fabs(eul*dIdt[genotype_idx][plasmid_idx]) > vanish_bound)
                {
                    converged = false;
                }

                I[genotype_idx][plasmid_idx] += eul*dSdt[genotype_idx];
            } // end for (int plasmid_idx
        } // end for (int genotype_idx = 0; genotype_idx < 4; ++genotype_idx)
   
        if (time_idx % skip_output == 0)
        {
            write_data(output_file);
        }

        if (converged)
        {
           break;
        } 
    } // end for int time_idx
} // end void iterateSys()
