#include "individual.hpp"

Individual::Individual() :
    p_chr{false}
    ,t_chr{false}
    ,p_plasmid{false}
    ,t_plasmid{false}
    ,has_plasmid{false}
{}

Individual::Individual(Individual const &other) :
    p_chr{other.p_chr}
    ,t_chr{other.t_chr}
    ,p_plasmid{other.p_plasmid}
    ,t_plasmid{other.t_plasmid}
    ,has_plasmid{other.has_plasmid}
{}

void Individual::operator=(Individual const &other)
{
    p_chr = other.p_chr;
    t_chr = other.t_chr;
    p_plasmid = other.p_plasmid;
    t_plasmid = other.t_plasmid;
    has_plasmid = other.has_plasmid;
}
