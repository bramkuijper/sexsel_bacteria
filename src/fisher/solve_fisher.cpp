#include <iostream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include "solve_fisher.hpp"


SolveFisher::SolveFisher(int argc, char **argv) 
    :
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
    ,max_time{0}
    ,gamma{0.0}
    ,d{0.0}
    ,pi{0.0}
    ,cp{0.0}
    ,ct{0.0}
    ,ht{0.0}
    ,hp{0.0}
    ,mu_t{0.0,0.0}
    ,mu_p{0.0,0.0}
    ,r{0.0}
    ,a{0.0}
    ,N{0.0}
    ,frac_infected_t0{0.0}
    ,p2_t0{0.0} 
    ,t2_t0{0.0} 
    ,D_t0{0.0} 
    ,skip_output{1000}
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
    gamma = atof(argv[2]);
    delta = atof(argv[3]);
    max_time = atof(argv[4]);
    d = atof(argv[5]);
    N = atof(argv[6]);
    frac_infected_t0 = atof(argv[7]);
    p2_t0 = atof(argv[8]);
    t2_t0 = atof(argv[9]);
    D_t0 = atof(argv[10]);
    cp = atof(argv[11]);
    ct = atof(argv[12]);
    pi = atof(argv[13]);
    r = atof(argv[14]);
    a = atof(argv[15]);
    hp = atof(argv[16]);
    ht = atof(argv[17]);
    mu_t[0] = atof(argv[18]);
    mu_t[1] = atof(argv[19]);
    mu_p[0] = atof(argv[20]);
    mu_p[1] = atof(argv[21]);

    base_name = argv[22];

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
        int nt2_recip_chr = geno == t2p1 || geno == t2p2;

        // number of p2s on recipient's chromosome (0 or 1)
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
            int nt2_donor_plm = geno_plasmid_donor == t2p1 
                || geno_plasmid_donor == t2p2;

            // count number of p2s on plasmid of donor
            int np2_donor_plm = geno_plasmid_donor == t1p2 
                || geno_plasmid_donor == t2p2;

            // fecundity of infected
            bI[geno_recip_chr_idx][geno_donor_plm_idx] = 
                exp(-cp*(np2_recip_chr+np2_donor_plm) 
                        -ct*(nt2_recip_chr+nt2_donor_plm) - delta);

            // now loop through chromosomal genotypes of infected
            // and calculate infection rates
            for (int geno_donor_chr_idx = 0; 
                    geno_donor_chr_idx < 4; ++geno_donor_chr_idx)
            {
                genotype geno_donor_chr = 
                    static_cast<genotype>(geno_donor_chr_idx);

                int nt2_donor_chr = geno_donor_chr == t2p1 || 
                    geno_donor_chr == t2p2;

                // by default set the infection rate to 0
                beta_SxI[geno_recip_chr_idx][
                    geno_donor_plm_idx][geno_donor_chr_idx] = 1.0;

                // if recipient is choosy on either chrom or plasmid
                // and something in donor that looks like an orament
                if (np2_recip_chr > 0 && nt2_donor_plm + nt2_donor_chr > 0)
                {
                    beta_SxI[geno_recip_chr_idx][geno_donor_plm_idx][geno_donor_chr_idx] += 
                        nt2_donor_plm + nt2_donor_chr == 1 ? 
                            ht * a // donor heterozygous t1 t2
                            :
                            a; // donor homozygous
                }

                for (int geno_recip_plm_idx = 0; 
                        geno_recip_plm_idx < 4; ++geno_recip_plm_idx)
                {
                    genotype geno_recip_plm = 
                        static_cast<genotype>(geno_recip_plm_idx);

                    int np2_recip_plm = geno_recip_plm == t1p2 || geno_recip_plm == t2p2;

                    // infection rate from the perspective of an infected recipient
                    // getting infected by an infected donor
                    beta_IxI[geno_recip_chr_idx][geno_recip_plm_idx][
                        geno_donor_chr_idx][geno_donor_plm_idx] = 1.0;

                    // recipient heterozygous for p2
                    // choosing any ornamented donor
                    if (np2_recip_chr + np2_recip_plm == 1 
                            && nt2_donor_plm + nt2_donor_chr > 0)
                    {
                        beta_IxI[geno_recip_chr_idx][geno_recip_plm_idx][
                            geno_donor_chr_idx][geno_donor_plm_idx] +=
                               nt2_donor_plm + nt2_donor_chr == 1 ?
                                ht * hp * a
                                :
                                hp * a;
                    }// recipient homozygous for p2
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


//    for (int idx1 = 0; idx1 < 4; ++idx1)
//    {
//        for (int idx2 = 0; idx2 < 4; ++idx2)
//        {
//            for (int idx3 = 0; idx3 < 4; ++idx3)
//            {
//                for (int idx4 = 0; idx4 < 4; ++idx4)
//                {
//                    std::cout << idx1 << " " << idx2 
//                        << " " << idx3 << " " << idx4 << " " << beta_IxI[idx1][idx2][idx3][idx4] << std::endl;
//                    if (idx1 < 2 && idx2 < 2)
//                    {
//                        assert(beta_IxI[idx1][idx2][idx3][idx4] == 1.0);
//                    }
//
//                }
//            }
//        }
//    }

    assert(beta_IxI[3][3][3][3] == 1 + a);

    // this would be a genotype t1p1 x t1p1 meeting any other genotype
    assert(beta_IxI[0][0][3][3] == 1.0);

}//end init_arguments

// write all the parameters to the file data_file
void SolveFisher::write_parameters(std::ofstream &data_file)
{
    for (int geno_chr_idx = 0; geno_chr_idx < 4; ++geno_chr_idx)
    {
        genotype geno_chr = static_cast<genotype>(geno_chr_idx);
        std::string tval_chr = geno_chr == t2p1 || geno_chr == t2p2 ? "t2" : "t1";
        std::string pval_chr = geno_chr == t1p2 || geno_chr == t2p2 ? "p2" : "p1";

        // write fecundity parameters for all suscetible genotypes
        data_file << "fec_S" 
            << tval_chr 
            << pval_chr 
            << ";" 
            << bS[geno_chr_idx]
            << std::endl;

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
                bI[geno_chr][geno_plm_idx]
                << std::endl;
        
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
        << "r;" << r << std::endl
        << "eul;" << eul << std::endl
        << "frac_infected_t0;" << frac_infected_t0 << std::endl
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

// write headers of the data file
void SolveFisher::write_data_headers(std::ofstream &data_file)
{
    data_file << "time;";

    for (int chr_idx = 0; chr_idx < 4; ++chr_idx)
    {
        data_file << "St" 
            << (1 + has_t2[chr_idx]) 
            << "p"
            << (1 + has_p2[chr_idx]) 
            << ";";

        for (int plm_idx = 0; plm_idx < 4; ++plm_idx)
        {
            data_file << "It"
                << (1 + has_t2[chr_idx]) 
                << "p"
                << (1 + has_p2[chr_idx]) 
                << "t"
                << (1 + has_t2[plm_idx]) 
                << "p"
                << (1 + has_p2[plm_idx]) 
                << ";";
        }
    }

    data_file << "S;I;N;"
                << "p2;p2_chr;p2_plm;"
                << "t2;t2_chr;t2_plm;";
    data_file << std::endl;
} // end SolveFisher::write_data_headers()

// writes data to output file
void SolveFisher::write_data(std::ofstream &data_file, int const time_step)
{
    data_file << time_step << ";";

    double total_S = 0.0;
    double total_I = 0.0;

    int p2 = 0;
    int p2_plm = 0;
    int p2_chr = 0;
    
    int t2 = 0;
    int t2_plm = 0;
    int t2_chr = 0;

    for (int chr_idx = 0; chr_idx < 4; ++chr_idx)
    {
        data_file << S[chr_idx] << ";";

        total_S += S[chr_idx];

        if (has_p2[chr_idx])
        {
            p2 += S[chr_idx];
            p2_chr += S[chr_idx];
        }

        if (has_t2[chr_idx])
        {
            t2 += S[chr_idx];
            t2_chr += S[chr_idx];
        }

        for (int plm_idx = 0; plm_idx < 4; ++plm_idx)
        {
            data_file << I[chr_idx][plm_idx] << ";";
            total_I += I[chr_idx][plm_idx];

            if (has_p2[plm_idx])
            {
                p2 += I[chr_idx][plm_idx];
                p2_plm += I[chr_idx][plm_idx];
            }

            if (has_p2[chr_idx])
            {
                p2 += I[chr_idx][plm_idx];
                p2_chr += I[chr_idx][plm_idx];
            }

            if (has_t2[plm_idx])
            {
                t2 += I[chr_idx][plm_idx];
                t2_plm += I[chr_idx][plm_idx];
            }

            if (has_t2[chr_idx])
            {
                t2 += I[chr_idx][plm_idx];
                t2_chr += I[chr_idx][plm_idx];
            }
        }
    }

    data_file << total_S << ";" 
        << total_I << ";" 
        << total_S + total_I << ";";

    data_file << (double)p2 / (total_S + 2*total_I) << ";";

    data_file << (double)p2_chr / (total_S + total_I) << ";";
    data_file << (double)p2_plm / total_I << ";";

    data_file << (double)t2 / (total_S + 2*total_I) << ";";
    data_file << (double)t2_chr / (total_S + total_I) << ";";
    data_file << (double)t2_plm / total_I << ";";

    data_file << std::endl;

} // SolveFisher::write_data()

// initialize the numbers in this population
void SolveFisher::init_population()
{
    // calculate genotype frequencies
    double freqs[4] = 
        {(1.0 - t2_t0) * (1.0 - p2_t0) - r * D_t0
            ,t2_t0 * (1.0 - p2_t0) + r * D_t0
            ,(1.0 - t2_t0) * p2_t0 + r * D_t0
            ,t2_t0 * p2_t0 - r * D_t0
        };

    for (int S_idx = 0; S_idx < 4; ++S_idx)
    {
        S[S_idx] = N * (1.0 - frac_infected_t0) * freqs[S_idx];
        
        assert(freqs[S_idx] >= 0);
        assert(freqs[S_idx] <= 1.0);
        assert(S[S_idx] >= 0);
        assert(S[S_idx] <= N * (1.0 - frac_infected_t0));
    }

    for (int I_idx_chr = 0; I_idx_chr < 4; ++I_idx_chr)
    {
        for (int I_idx_plm = 0; I_idx_plm < 4; ++I_idx_plm)
        {
            I[I_idx_chr][I_idx_plm] = 
                N * frac_infected_t0 * freqs[I_idx_chr] * freqs[I_idx_plm];

            assert(I[I_idx_chr][I_idx_plm] >= 0);
            assert(I[I_idx_chr][I_idx_plm] <= N * frac_infected_t0);
        }
    }
}// end SolveFisher::init_population()

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

    write_data_headers(output_file);

    init_population();

    int Ntplus1;

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
        double recombination_in[4][4];
        double recombination_out[4][4];

        double total_recombination = 0.0;

        // - calculate plasmid loss rate
        // - reset force of infection calculation
        for (int genotype_chr_idx = 0; genotype_chr_idx < 4; ++genotype_chr_idx)
        {
            assert(S[genotype_chr_idx] >= 0);
            assert(S[genotype_chr_idx] <= N);

            total_force_of_infection_lossS[genotype_chr_idx] = 0.0;

            // level 2: donor plasmid 
            for (int genotype_plm_idx = 0; 
                    genotype_plm_idx < 4; ++genotype_plm_idx)
            {
                assert(I[genotype_chr_idx][genotype_plm_idx] >= 0);
                assert(I[genotype_chr_idx][genotype_plm_idx] <= N);

                // set total rate of gain in I due to infections to 0 and
                // calculate value in the next loop 
                total_force_of_infection_gainI[genotype_chr_idx][genotype_plm_idx] = 0.0;

                sum_plasmid_loss[genotype_chr_idx] += 
                    gamma * I[genotype_chr_idx][genotype_plm_idx];
            
                recombination_in[genotype_chr_idx][genotype_plm_idx] = 
                    recombination_out[genotype_chr_idx][genotype_plm_idx] = 0.0;
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

                    // influx of infecteds with genotype 
                    // I[geno_recip_chr_idx][geno_recip_plm_idx]
                    // due to recombination
                    //
                    // note the 1/2 there, as only half of the recombinants will 
                    // have the desired combination of [tx py] [tz pv] on chromosome
                    // and plasmid respectively, the other half will have 
                    // [tz pv] [tx py] on chromosome and plasmid respectively
                    recombination_in[geno_recip_chr_idx][geno_recip_plm_idx] +=
                        0.5 * r * (I[recombinant1][recombinant2] +
                                I[recombinant2][recombinant1]);

                    recombination_out[geno_recip_chr_idx][geno_recip_plm_idx] += 
                        r * I[geno_recip_chr_idx][geno_recip_plm_idx];

                    total_recombination += 
                        recombination_in[geno_recip_chr_idx][geno_recip_plm_idx] 
                        - recombination_out[geno_recip_chr_idx][geno_recip_plm_idx];
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

        if (total_recombination >= 1e-07)
        {
            std::cout << time_idx << " " << total_recombination << std::endl;
        }


        assert(fabs(total_recombination) < 1e-07);

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
            total_mutations += mu_t[!is_t2] * S[geno_chr_oppT] + mu_p[!is_p2] * S[geno_chr_oppP];
            total_mutations -= (mu_t[is_t2] + mu_p[is_p2]) * S[genotype_idx];

            // dSdt
            dSdt[genotype_idx] = 
                // 1. birth of new susceptibles
                bS[genotype_idx] * (1.0 - kappa * N) * S[genotype_idx] +

                // 2. gain in susceptibles due to loss of plasmid
                sum_plasmid_loss[genotype_idx] +

                // 3. mutation incoming
                + mu_t[!is_t2] * S[geno_chr_oppT] + mu_p[!is_p2] * S[geno_chr_oppP]

                // 4. mutation outgoing
                - (mu_t[is_t2] + mu_p[is_p2]) * S[genotype_idx] 
                
                // 5. loss due to deaths
                - d * S[genotype_idx] 
                
                // 6. loss due to infections
                - (1.0-pi) * total_force_of_infection_lossS[genotype_idx] * S[genotype_idx]; 

            for (int plasmid_idx = 0; plasmid_idx < 4; ++plasmid_idx)
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
                    + (1.0-pi) * total_force_of_infection_gainI[genotype_idx][plasmid_idx] *
                        I[genotype_idx][plasmid_idx]

                    // 4. loss in infecteds
                    - (d + gamma) * I[genotype_idx][plasmid_idx]

                    // 5. recombination
                    + recombination_in[genotype_idx][plasmid_idx]
                    - recombination_out[genotype_idx][plasmid_idx];

            } // end for plasmid_idx

        } // end for int genotype_idx

        converged = true;

        N = 0;
            
//        if (time_idx > 20000)
//        {
//            std::cout << "t;" << time_idx << ";";
//        }

        for (int genotype_idx = 0; genotype_idx < 4; ++genotype_idx)
        {
            // S[i](t+1) = S[i] + eul * dSdt[i]
            // S[i](t+1) - S[i] < vanish_bound
            // eul*dSdt < vanish_bound
            if (fabs(eul*dSdt[genotype_idx]) > vanish_bound)
            {
                converged = false;
            }

//            if (time_idx > 20000)
//            {
//                std::cout << "S" << genotype_idx << ";" << eul * dSdt[genotype_idx] << ";";
//            }

            // update value of S[i]
            S[genotype_idx] = S[genotype_idx] + eul*dSdt[genotype_idx];

            if (S[genotype_idx] < 0.0)
            {
                S[genotype_idx] = 0.0;
            }

            N += S[genotype_idx];

            for (int plasmid_idx = 0; plasmid_idx < 4; ++plasmid_idx)
            {
                if (fabs(eul*dIdt[genotype_idx][plasmid_idx]) > vanish_bound)
                {
                    converged = false;
                }

                I[genotype_idx][plasmid_idx] += eul*dIdt[genotype_idx][plasmid_idx];

                if (I[genotype_idx][plasmid_idx] < 0.0)
                {
                    I[genotype_idx][plasmid_idx] = 0.0;
                }
                
//                if (time_idx > 20000)
//                {
//                    std::cout << "I" << genotype_idx << "_" << plasmid_idx << ";" << eul * dIdt[genotype_idx][plasmid_idx] << ";";
//                }

                N += I[genotype_idx][plasmid_idx];
            } // end for (int plasmid_idx

        } // end for (int genotype_idx = 0; genotype_idx < 4; ++genotype_idx)
   
//        if (time_idx > 20000)
//        {
//            std::cout << std::endl;
//        }
        if (time_idx % skip_output == 0)
        {
            write_data(output_file, time_idx);
        }

        if (converged)
        {
           break;
        } 


    } // end for int time_idx
} // end void iterateSys()
