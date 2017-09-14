#pragma once

//#include "Global.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "Vec2d.h"
#include "Particle.h"

using namespace std;

static inline bool sort_pid(const Particle &i, const Particle &j) {return i.pid < j.pid;}

class NonAffine {
 private:
  int Nx, Ny, N;
  const int num_nebz = 6;
  Vec2d l;
  
  std::vector<Particle> ref_pos;
  std::vector<int> nebz_pid;
  std::vector<std::vector<int>> nebz2_pid;

  double mean_chi;
  double hchi;

  std::vector<double> chi_buff;

  // X_prev and eps_prev of the particle to be displaced
  vector<double> X_prev;
  vector<double> eps_prev;
  vector<int> get_nebz1(const int i);

 public:
  
   
  void calc_X(vector<Particle> &particles);

  void calc_Y(vector<Particle> &particles);
  void calc_eps(vector<Particle> &particles);

  // save & restore
  void save_X_eps(int i);
  void restore_X_eps(int i);

  // Move these to private
  std::vector<double> Y;
  std::vector<double> X;
  std::vector<double> eps;

  void init(int,int,Vec2d,double);
  void set_ref(const std::vector<Particle>&);
  void scale_ref(const Vec2d& s);

  double calc_chi(vector<Particle>& particles);



  // strain ref
  void apply_strain_on_ref_lat(const double e) {
    for (size_t i = 0; i < N; i++) {
      ref_pos[i].r.x = (1.0 + e) * ref_pos[i].r.x;
      ref_pos[i].r.y = (1.0 - e) * ref_pos[i].r.y;
    }
  }
  
  // PBC
	void pbc(Vec2d &r) {
		r.x -= l.x * rint(r.x / l.x);
		r.y -= l.y * rint(r.y / l.y);
	}
	
	inline double pbcx(double xx) { return (xx - l.x * rint(xx/l.x) );}
	inline double pbcy(double yy) { return (yy - l.y * rint(yy/l.y) );}
};
