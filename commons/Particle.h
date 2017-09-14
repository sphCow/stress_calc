#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Vec2d.h"

struct Particle{
    int pid;
    Vec2d r;
		Vec2d r_pbc;
		double chi;
        double chi2;
    double pe;
		double virial_tensor[4];
    int cid;
    Vec2d f;
    
    double virial_na[4];
    double virial_na2[4];
    double virial_nb[4];
    
};

struct Particle_Ising{
    int pid;
    int spin;
};

struct Header {
	double bin_id;
	Vec2d l;
};

#endif
