#ifndef PLASMID_HPP
#define PLASMID_HPP


class Plasmid 
{
    public:
    
        // genotype	
        //preference in plasmid 
        //trait in plasmid
        bool p_gen;
        bool t_gen;

        Plasmid();

        Plasmid(Plasmid const &other);

        void operator=(Plasmid const &other);


};

#endif
