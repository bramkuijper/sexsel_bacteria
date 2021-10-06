#include "plasmid.hpp"

Plasmid::Plasmid() :
    ,p_gen{false}
    ,t_gen{false}
{}

Plasmid::Plasmid(Plasmid const &other) :
    ,p_gen{other.p_gen}
    ,t_gen{other.t_gen}
{}

void Plasmid::operator=(Plasmid const &other)
{
    p_gen = other.p_gen;
    t_gen = other.t_gen;
}
