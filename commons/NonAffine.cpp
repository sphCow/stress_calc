#include "NonAffine.h"

vector<int> NonAffine::get_nebz1(const int i) {
  vector<int> v;
  v.resize(6);
  
  int ix = int(i%Nx);
  int iy = int(i/Ny);

  if (iy%2 == 0) {
    v[0] = ((ix + Nx - 1) % Nx + ((iy + Ny + 1) % Ny) * Nx);  // NW
    v[3] = ((ix + Nx - 1) % Nx + ((iy + Ny - 1) % Ny) * Nx);  // SW
  }
  
  else {
    v[0] = ((ix + 1) % Nx + ((iy + Ny + 1) % Ny) * Nx);  // NE
    v[3] = ((ix + 1) % Nx + ((iy + Ny - 1) % Ny) * Nx);  // SE
  }
  
  v[2] = ((ix + 1) % Nx + ((iy + Ny) % Ny) * Nx);  // EAST
  v[1] = ((ix + Nx) % Nx + ((iy + Ny + 1) % Ny) * Nx);  // N
  v[4] = ((ix + Nx - 1) % Nx + ((iy + Ny) % Ny) * Nx);  // W
  v[5] = ((ix + Nx) % Nx + ((iy + Ny - 1) % Ny) * Nx);  // S

  return v;
}

void NonAffine::init(int Nx_, int Ny_, Vec2d l_, double h_) {
  // ref_pos.resize(N);
  Nx = Nx_;
  Ny = Ny_;
  l = l_;
  hchi = h_;
  N = Nx*Ny;
  
  Y.resize(4 * N);
  X.resize(4 * N);
  eps.resize(4 * N);
  nebz_pid.resize(num_nebz*N);
  nebz2_pid.resize(N);

  chi_buff.resize(num_nebz + 1);
  X_prev.resize(4 * (num_nebz + 1));
  eps_prev.resize(4 * (num_nebz + 1));
  
  
  
  
}

void NonAffine::set_ref(const std::vector<Particle>& _particles) {
  ref_pos = _particles;
  
  // populate nebz
  for (long ix = 0; ix < Nx; ix++) {
    for (long iy = 0; iy < Ny; iy++) {
      long i = ix + iy * Nx;
      // clang-format off
      if (iy % 2 == 0) {
        nebz_pid[num_nebz * i + 0] = ((ix + Nx - 1) % Nx + ((iy + Ny + 1) % Ny) * Nx);  // NW
        nebz_pid[num_nebz * i + 3] = ((ix + Nx - 1) % Nx + ((iy + Ny - 1) % Ny) * Nx);  // SW
      }
      
      else {
        nebz_pid[num_nebz * i + 0] = ((ix + 1) % Nx + ((iy + Ny + 1) % Ny) * Nx);  // NE
        nebz_pid[num_nebz * i + 3] = ((ix + 1) % Nx + ((iy + Ny - 1) % Ny) * Nx);  // SE
      }
      
      nebz_pid[num_nebz * i + 2] = ((ix + 1) % Nx + ((iy + Ny) % Ny) * Nx);  // EAST
      nebz_pid[num_nebz * i + 1] = ((ix + Nx) % Nx + ((iy + Ny + 1) % Ny) * Nx);  // N
      nebz_pid[num_nebz * i + 4] = ((ix + Nx - 1) % Nx + ((iy + Ny) % Ny) * Nx);  // W
      nebz_pid[num_nebz * i + 5] = ((ix + Nx) % Nx + ((iy + Ny - 1) % Ny) * Nx);  // S
      // clang-format on
    }
  }
  
  // 2nd nebz without the 1st nebz
  for(int i=0; i<N; i++) {
    vector<int> nn1 = get_nebz1(i);
    
    for(auto &k:nn1) {
      vector<int> nn2 = get_nebz1(k);
      
      for(auto &m:nn2) {
        nebz2_pid[i].push_back(m);
      }
    }
  }
  
  //erase first neb
  for(int i=0; i<N; i++) {
    vector<int> nn1 = get_nebz1(i);
    
    for(auto &k:nn1) {
      nebz2_pid[i].erase(remove(nebz2_pid[i].begin(), nebz2_pid[i].end(), k), nebz2_pid[i].end());
    }
    
    //unique
    sort(nebz2_pid[i].begin(), nebz2_pid[i].end()); // 1 1 2 2 3 3 3 4 4 5 5 6 7 
    auto last = std::unique(nebz2_pid[i].begin(), nebz2_pid[i].end());
    // v now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
    nebz2_pid[i].erase(last, nebz2_pid[i].end());
    
  }
  
  //cout << "2nd nebz populated" << endl;
  //ofstream ofn;
  //ofn.open("nebz.txt");
  //for(int i=0; i<N; i++) {
  //  ofn << i << " -1 " << ref_pos[i].r.x << " " << ref_pos[i].r.y << endl;
  //  for(int j=0; j<nebz2_pid[i].size(); j++) {
  //    ofn << j << " " << nebz2_pid[i][j] << " " << ref_pos[nebz2_pid[i][j]].r.x << " " << ref_pos[nebz2_pid[i][j]].r.y << endl;
  //  }
  //  ofn << endl << endl;
  //}
  
  calc_Y(ref_pos);
}

void NonAffine::scale_ref(const Vec2d& s) {
  for (auto& i : ref_pos) {
    i.r.x += s.x * i.r.x;
    i.r.y += s.y * i.r.y;
  }
}

// Do this after setting ref lattice
void NonAffine::calc_Y(vector<Particle> &particles) {
  for (auto& pri : particles) {
    // pri -> i-th particle in current config
    // Y[4*i+j] gives me j-th component of Y matrix for i-th particle

    int i = pri.pid;
    Y[4 * i + 0] = 0.0;
    Y[4 * i + 1] = 0.0;
    Y[4 * i + 2] = 0.0;
    Y[4 * i + 3] = 0.0;

    Particle pRi = ref_pos[i];  // pRi -> i-th particle in reference lattice

    for (int j = 0; j < num_nebz; j++) {  // LOOP OVER NEBZ OF i

      int k = nebz_pid[num_nebz * i + j];

      Particle pRj = ref_pos[k];  // pRj -> j-th neighbour of i the particle in
                                  // reference lattice
      Particle prj = particles[k];  // prj -> j-th neighbour of i the particle
                                    // in current config

      Y[4 * i + 0] += pbcx(pRj.r.x - pRi.r.x) * pbcx(pRj.r.x - pRi.r.x);
      Y[4 * i + 1] += pbcx(pRj.r.x - pRi.r.x) * pbcy(pRj.r.y - pRi.r.y);
      Y[4 * i + 2] += pbcy(pRj.r.y - pRi.r.y) * pbcx(pRj.r.x - pRi.r.x);
      Y[4 * i + 3] += pbcy(pRj.r.y - pRi.r.y) * pbcy(pRj.r.y - pRi.r.y);
    }
  }
}

void NonAffine::calc_X(vector<Particle> &particles) {
  for (auto& pri : particles) {  // LOOP OVER ALL PARTICLES
    int i = pri.pid;
    X[4 * i + 0] = 0.0;
    X[4 * i + 1] = 0.0;
    X[4 * i + 2] = 0.0;
    X[4 * i + 3] = 0.0;

    // pri -> i-th particle in current config
    Particle pRi = ref_pos[i];  // pRi -> i-th particle in reference lattice

    // CALCULATE X
    for (int j = 0; j < num_nebz; j++) {  // num_nebz = 6 is hard coded
      int k = nebz_pid[num_nebz * i + j];
      Particle pRj = ref_pos[k];
      Particle prj = particles[k];

      X[4 * i + 0] += pbcx(prj.r.x - pri.r.x) * pbcx(pRj.r.x - pRi.r.x);
      X[4 * i + 1] += pbcx(prj.r.x - pri.r.x) * pbcy(pRj.r.y - pRi.r.y);
      X[4 * i + 2] += pbcy(prj.r.y - pri.r.y) * pbcx(pRj.r.x - pRi.r.x);
      X[4 * i + 3] += pbcy(prj.r.y - pri.r.y) * pbcy(pRj.r.y - pRi.r.y);
    }
  }
}

void NonAffine::calc_eps(vector<Particle> &particles) {
  for (auto& pri : particles) {  // LOOP OVER ALL PARTICLES
    int i = pri.pid;

    // epsilon = XY^-1 - I
    double inv_det_Y = 1.0 / (-Y[4 * i + 1] * Y[4 * i + 2] + Y[4 * i + 0] * Y[4 * i + 3]);  // |Y|

    eps[4 * i + 0] = inv_det_Y * (-X[4 * i + 1] * Y[4 * i + 2] +X[4 * i + 0] * Y[4 * i + 3]) -1.0;
    eps[4 * i + 1] = inv_det_Y * (X[4 * i + 1] * Y[4 * i + 0] - X[4 * i + 0] * Y[4 * i + 1]);
    eps[4 * i + 2] = inv_det_Y * (-X[4 * i + 3] * Y[4 * i + 2] + X[4 * i + 2] * Y[4 * i + 3]);
    eps[4 * i + 3] = inv_det_Y * (X[4 * i + 3] * Y[4 * i + 0] - X[4 * i + 2] * Y[4 * i + 1]) -1.0;
  }
}

/**********************************************************************/
/************************* Calculate Chi ******************************/
/**********************************************************************/

/*********************** CALC CHI FOR ALL PARTICLES *******************/
/* SORT PID was removed*/
double NonAffine::calc_chi(vector<Particle> &particles) {
  sort(particles.begin(), particles.end(), sort_pid);
  calc_X(particles);
  calc_eps(particles);

  mean_chi = 0.0;

  for (auto& pri : particles) {  // LOOP OVER ALL PARTICLES
    vector<int> traversed_nebz(36);
    
    double chi = 0.0;
    for(int j=0; j<4; j++) pri.virial_na[j] = 0.0;
    for(int j=0; j<4; j++) pri.virial_na2[j] = 0.0;

    int i = pri.pid;
    Particle pRi = ref_pos[i];  // pRi -> i-th particle in reference lattice

    // CHI
    Vec2d fhchi = Vec2d(0, 0);  // initilized to zero by default, i.e. hf.x = 0; hf.y = 0
    for (int j = 0; j < num_nebz; j++) {
      int k = nebz_pid[num_nebz * i + j];  // this is the index of the j-th neighbour of i-th particle

      Particle pRj = ref_pos[k];
      Particle prj = particles[k];

      // Vec2d dr(pbcx(prj.r.x - pri.r.x), pbcy(prj.r.y - pri.r.y));
      // Vec2d dR(pbcx(pRj.r.x - pRi.r.x), pbcy(pRj.r.y - pRi.r.y));

      Vec2d dr, dR;
      dr.x = pbcx(pri.r.x - prj.r.x);
      dr.y = pbcy(pri.r.y - prj.r.y);

      dR.x = pbcx(pRi.r.x - pRj.r.x);
      dR.y = pbcy(pRi.r.y - pRj.r.y);

      double chi_x = dr.x - (1.0 + eps[4 * i + 0]) * dR.x - (eps[4 * i + 1]) * dR.y;
      double chi_y = dr.y - (eps[4 * i + 2]) * dR.x - (1.0 + eps[4 * i + 3]) * dR.y;
      chi += pow(chi_x, 2.0) + pow(chi_y, 2.0);

      fhchi.x += 4.0*dr.x - 2.0*(1.0 + eps[4*k+0])*dR.x - 2.0*(eps[4*k+1])*dR.y;
      fhchi.y += 4.0*dr.y - 2.0*(      eps[4*k+2])*dR.x - 2.0*(1.0 + eps[4*k+3])*dR.y;
      
      //pri.virial_na[0] += (2*dr.x - dR.x - eps[4*k+0]*dR.x - eps[4*k+1]*dR.y)*dr.x;
      //pri.virial_na[1] += (2*dr.y - dR.y - eps[4*k+2]*dR.x - eps[4*k+3]*dR.y)*dr.y;
      
    }

    pri.chi2 = chi;    // chi_pp array stores the per particle values of chi
    mean_chi += chi;  // accumulates chi_i

    pri.f += (fhchi)*(hchi); // add
    //pri.virial_na[0] *= (-2*hchi);
    //pri.virial_na[1] *= (-2*hchi);
    
    pri.virial_na[0] = pri.f.x*pri.r.x;
    pri.virial_na[1] = pri.f.y*pri.r.y;
    
  }
  
  // comapre with SUSE's expression
  for(int i=0; i<N; i++) {
    Particle pri = particles[i];
    
    //first term
    for (int j = 0; j < num_nebz; j++) {
      int k = nebz_pid[num_nebz * i + j];  // this is the index of the j-th neighbour of i-th particle
      
      Particle prj = particles[k];
      
      Vec2d dr;
      dr.x = pbcx(pri.r.x - prj.r.x);
      dr.y = pbcy(pri.r.y - prj.r.y);
      
      
      particles[i].virial_na2[0] += 7*dr.x*dr.x;
      particles[i].virial_na2[1] += 7*dr.y*dr.y;
    }
    
    //first term
    for (int j = 0; j < nebz2_pid[i].size(); j++) {
      int k = nebz2_pid[i][j];  // this is the index of the j-th neighbour of i-th particle
      
      Particle prj = particles[k];
      
      Vec2d dr;
      dr.x = pbcx(pri.r.x - prj.r.x);
      dr.y = pbcy(pri.r.y - prj.r.y);
      
      
      particles[i].virial_na2[0] -= dr.x*dr.x;
      particles[i].virial_na2[1] -= dr.y*dr.y;
    }
    
    // put in factors
    particles[i].virial_na2[0] = particles[i].virial_na2[0]*hchi/6.0;
    particles[i].virial_na2[1] = particles[i].virial_na2[1]*hchi/6.0;
    
  }

  mean_chi /= double(N);
  return mean_chi;
}



/**********************************************************************/
/************************* for i-th particle **************************/
/**********************************************************************/

/* calc_X & calc_eps functions updates X & eps matrix of i-th particle
 * and its neighbours and saves the previous values of X & eps mat
 * of i and upto its 2nd nn */

/*
void NonAffine::calc_X(Particle& pri) {
  int i = pri.pid;

  // save X[i] and X[{nn of i}]
  X_prev[0] = X[4 * i];
  X_prev[1] = X[4 * i + 1];
  X_prev[2] = X[4 * i + 2];
  X_prev[3] = X[4 * i + 3];

  for (int j = 0; j < num_nebz; j++) {
    int neb = nebz_pid[num_nebz * i + j];

    for (int k = 0; k < 4; k++) X_prev[4 + 4 * j + k] = X[4 * neb + k];
  }

  X[4 * i + 0] = 0.0;
  X[4 * i + 1] = 0.0;
  X[4 * i + 2] = 0.0;
  X[4 * i + 3] = 0.0;

  // pri -> i-th particle in current config
  Particle pRi = ref_pos[i];  // pRi -> i-th particle in reference lattice

  // CALCULATE X
  for (int j = 0; j < num_nebz; j++) {  // num_nebz = 6 is hard coded
    int k = nebz_pid[num_nebz * i + j];
    Particle pRj = ref_pos[k];
    Particle prj = particles[k];

    X[4 * i + 0] += pbcx(prj.r.x - pri.r.x) * pbcx(pRj.r.x - pRi.r.x);
    X[4 * i + 1] += pbcx(prj.r.x - pri.r.x) * pbcy(pRj.r.y - pRi.r.y);
    X[4 * i + 2] += pbcy(prj.r.y - pri.r.y) * pbcx(pRj.r.x - pRi.r.x);
    X[4 * i + 3] += pbcy(prj.r.y - pri.r.y) * pbcy(pRj.r.y - pRi.r.y);
  }
}

void NonAffine::calc_eps(Particle& pri) {
  int i = pri.pid;

  // save X[i] and X[{nn of i}]
  eps_prev[0] = eps[4 * i];
  eps_prev[1] = eps[4 * i + 1];
  eps_prev[2] = eps[4 * i + 2];
  eps_prev[3] = eps[4 * i + 3];

  for (int j = 0; j < num_nebz; j++) {
    int neb = nebz_pid[num_nebz * i + j];

    for (int k = 0; k < 4; k++) eps_prev[4 + 4 * j + k] = eps[4 * neb + k];
  }

  // epsilon = XY^-1 - I
  // clang-format off
  double inv_det_Y = 1.0/(-Y[4*i+1]*Y[4*i+2]+Y[4*i+0]*Y[4*i+3]);  // |Y|

  eps[4*i+0] = inv_det_Y*(-X[4*i+1]*Y[4*i+2]+X[4*i+0]*Y[4*i+3])-1.0;
  eps[4*i+1] = inv_det_Y*(X[4*i+1]*Y[4*i+0]-X[4*i+0]*Y[4*i+1]);
  eps[4*i+2] = inv_det_Y*(-X[4*i+3]*Y[4*i+2]+X[4*i+2]*Y[4*i+3]);
  eps[4*i+3] = inv_det_Y*(X[4*i+3]*Y[4*i+0]-X[4*i+2]*Y[4*i+1]) -1.0;
  // clang-format on
}


double NonAffine::calc_chi(Particle& pri) {
  // sort(particles.begin(), particles.end(), sort_pid);
  calc_X(pri);
  calc_eps(pri);

  double chi = 0.0;

  int i = pri.pid;
  Particle pRi = ref_pos[i];  // pRi -> i-th particle in reference lattice

  // CHI
  Vec2d fhchi = Vec2d(0, 0);

  for (int j = 0; j < num_nebz; j++) {
    int k = nebz_pid[num_nebz * i + j];  // this is the index of the j-th
                                         // neighbour of i-th particle

    Particle pRj = ref_pos[k];
    Particle prj = particles[k];

    Vec2d dr(pbcx(pri.r.x - prj.r.x), pbcy(pri.r.y - prj.r.y));
    Vec2d dR(pbcx(pRi.r.x - pRj.r.x), pbcy(pRi.r.y - pRj.r.y));

    double chi_x =
        dr.x - (1.0 + eps[4 * i + 0]) * dR.x - (eps[4 * i + 1]) * dR.y;
    double chi_y =
        dr.y - (eps[4 * i + 2]) * dR.x - (1.0 + eps[4 * i + 3]) * dR.y;
    chi += pow(chi_x, 2.0) + pow(chi_y, 2.0);
  }

  pri.chi = chi;  // chi_pp array stores the per particle values of chi
  return chi;
}

// Calc Chi Partial - calc chi of particle id=rpos & its nn 
double NonAffine::calc_chi_partial(int rpos) {
  double new_nf_sum = calc_chi(particles[rpos]);  // chi for particle at rpos
  for (int j = 0; j < 6; j++) {
    int pid = particles[rpos].pid;
    int neb = nebz_pid[6 * pid + j];
    new_nf_sum += calc_chi(particles[neb]);
  }
  return new_nf_sum;
}

//********************* Restore *****************
void NonAffine::restore_X_eps(int i) {
  // restore X[i] and X[{nn of i}]
  for (int k = 0; k < 4; k++) {
    eps[4 * i + k] = eps_prev[k];
    X[4 * i + k] = X_prev[k];
  }

  for (int j = 0; j < num_nebz; j++) {
    int neb = nebz_pid[num_nebz * i + j];

    for (int k = 0; k < 4; k++) {
      eps[4 * neb + k] = eps_prev[4 + 4 * j + k];
      X[4 * neb + k] = X_prev[4 + 4 * j + k];
    }
  }
}
*/
