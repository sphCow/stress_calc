#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Particle.h"

using namespace std;

#ifndef NONBONDEDFORCESSTRESS
#define NONBONDEDFORCESSTRESS

class NonBondedForces {
private:
	int N;
	Vec2d l;
	Vec2d lh;
	Vec2d lcell;
	int ncell_x, ncell_y, ncell;
	vector<int> clen;
	vector<int> cbegin;
	void update_cell_list();
	static inline bool sort_cid(const Particle &i, const Particle &j) {return i.cid < j.cid;}
	static inline bool sort_pid(const Particle &i, const Particle &j) {return i.pid < j.pid;}
  void get_pair_LJ(Particle &i, Particle &j);
	double rc,rc2, c0, c2, c6;
	
	double pe;
	
	// PBC
	void pbc(Vec2d &r) {
		r.x -= l.x * rint(r.x / l.x);
		r.y -= l.y * rint(r.y / l.y);
	}
	
	inline double pbcx(double xx) { return (xx - l.x * rint(xx/l.x) );}
	inline double pbcy(double yy) { return (yy - l.y * rint(yy/l.y) );}
	
	vector<Particle> particles;
			
public:
	void init(double rc_, Vec2d l_, int N_, vector<Particle>& particles_);
	void calc();
	
	void get_results(vector<Particle> &particles_) {
		sort(particles.begin(), particles.end(), sort_pid);
		
		for(auto &ip: particles) {
			ip.virial_tensor[0] *= 0.5; 
			ip.virial_tensor[1] *= 0.5;
			ip.virial_tensor[2] *= 0.5;
			ip.virial_tensor[3] *= 0.5;
		}
		
		particles_ = particles;
		
	}
	
	double get_total_pe() {
		return 0.5*pe;
	}
	
};

#endif
