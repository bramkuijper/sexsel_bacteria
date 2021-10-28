#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP


class Individual 
{
    public:
    
        // genotype	
        //preference in chromosome 
        //trait in chromosome
        bool p_chr;
        bool t_chr;

        //preference in plasmid 
        //trait in plasmid
        bool p_plasmid;
        bool t_plasmid;
	
        // plasmid absent/present
        bool has_plasmid;

        Individual();

        Individual(Individual const &other);

        void operator=(Individual const &other);


};

#endif
