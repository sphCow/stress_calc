#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdio>
#include "Particle.h"
#include "nonbonded_stress.h"
#include "NonAffine.h"

#include "rapidjson/document.h"

using namespace rapidjson;
using namespace std;

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
  if(from.empty())
    return;
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

int main(int argc, char** argv) {
  string path = argv[1];
  
  /************ PARSE JSON ****************/
  string jsonstr,line;
  
  // Read JSON
  ifstream f;
  string jsonpath = path+"/c000000";
  f.open(jsonpath.c_str());
  
  while (getline(f, line)) {
    if(line.find("## INPUT JSON ##") != string::npos)
      continue;
    if(line.find("## END OF INPUT JSON ##") != string::npos)
      break;
    
    jsonstr += line;
    jsonstr += '\n';
  }
  
  f.close();
  //cout << "Parsing JSON block..." << endl;
  
  //remove #
  replaceAll(jsonstr, "#", "");
  
  Document d;
  d.Parse<kParseCommentsFlag>(jsonstr.c_str());
  
  // parse require
  int Nx = d["N"][0].GetInt();
  int Ny = d["N"][1].GetInt();
  int N = Nx*Ny;
  double rc = d["nb_rc"].GetDouble();
  double rho = d["rho"].GetDouble();
  double eps = d["epsilon"].GetDouble();
  double hchi = d["hchi"].GetDouble();
  Vec2d l,lh,a;
  
  //cout << N << " " << rc << " " << hchi << endl;
  
  //*********************** Generate reference ***********************//
  vector<Particle> pref;
  pref.resize(N);
  
  a.x = sqrt(2.0 / (rho * sqrt(3.0)));
  a.y = a.x * sqrt(3.0) / 2.0;
  
  l.x = Nx * a.x;
  l.y = Ny * a.y;
  lh = l * 0.5;
  
  // interchanged i & j [EDIT]
  for (long i = 0; i < Nx; i++) {
    for (long j = 0; j < Ny; j++) {
      double xx = -lh.x + i * a.x;
      double yy = -lh.y + j * a.y;
      
      if (j % 2 != 0) xx += a.x / 2.0;
      
      long id = i + j * Nx;
      
      pref[id].pid = id;
      pref[id].r.x = xx;  // x
      pref[id].r.y = yy;  // y
    }
  }
  
  //cout << "eps -> " << eps << endl;
  
  l.x += eps*l.x ;
  l.y -= eps*l.y ;
  
  for(int i=0; i<N; i++) {
    pref[i].r.x += eps*pref[i].r.x; // x
    pref[i].r.y -= eps*pref[i].r.y; 
  }
  
  string refpath = path+"ref";
  ofstream fref(refpath);
  for(int i=0; i<N; i++) {
    fref << pref[i].pid << " " << pref[i].r.x << " " << pref[i].r.y << endl;
  }
  
  /********************** Init NonAffine ******************/
  NonAffine nf;
  nf.init(Nx,Ny,l,hchi);
  nf.set_ref(pref);
  
  vector<Particle> particles;
  particles.resize(N);
  
  NonBondedForces nb;
  
  /*********************** Read Particles ******************/
  string header_str,ignore;
  
  int confid_int=stoi(argv[2]);
  char confid[64];
  snprintf(confid,64,"/c%06d",confid_int);
  
  string confpath = path+string(confid);
  f.open(confpath);
  int bin_id,pcount=0;
  double X_bin, bin_id_double;
  
  getline(f, line);
  header_str = line;
  
  //if(line.find("INTERMIDIATE") != std::string::npos) {
  //  cout << line << " -> " << "int bin" << endl;
  //  return EXIT_SUCCESS;
  //}
  
  replaceAll(line, "INTERMIDIATE", "");
  replaceAll(line, "#", "");
  replaceAll(line, "NEW_BIN", "");
  
  stringstream iss(line);
  iss >> bin_id_double >> ignore >> X_bin >> ignore >> l.x >> l.y;
  bin_id = bin_id_double*2;
  
  //cout << "header -> " << bin_id << " " << X_bin << " " << l.x << " " << l.y << endl;
  
  while (getline(f, line)) {
    stringstream iss(line);
    iss >> particles[pcount].pid 
        >> particles[pcount].r.x
        >> particles[pcount].r.y
        >> particles[pcount].chi;
    
    pcount++;
      
  }
    
  f.close();
  
  double X = nf.calc_chi(particles);
  
  nb.init(rc, l, N, particles);
  nb.calc();
  nb.get_results(particles);
  
  
  ofstream of(confpath+"_wstress", ios::out);
  of << header_str << endl;
  for(auto &ip:particles) {
    of << ip.pid << " " 
    << ip.r.x << " " 
    << ip.r.y << " " 
    << ip.chi << "\t"  
    << ip.virial_tensor[0]+ip.virial_na2[0] << " " << ip.virial_tensor[3]+ip.virial_na2[1] << " "
    << ip.virial_na2[0] << " " << ip.virial_na2[1] << endl;
  }
  
  of.close();
  
  // global stress 
  Vec2d sigma = Vec2d(0,0);
  //Vec2d sigma_na2 = Vec2d(0,0);
  //Vec2d sigma_nb = Vec2d(0,0);
  
  for(auto &ip:particles) {
    sigma.x += ip.virial_na2[0]+ip.virial_tensor[0];
    sigma.y += ip.virial_na2[1]+ip.virial_tensor[3];
    //sigma_na2.x += ip.virial_na2[0];
    //sigma_na2.y += ip.virial_na2[1];
    //sigma_nb.x += ip.virial_tensor[0];
    //sigma_nb.y += ip.virial_tensor[3];
  }
  
  sigma.x /= (l.x*l.y);
  sigma.y /= (l.x*l.y);
  //sigma_nb.x /= (l.x*l.y);
  //sigma_nb.y /= (l.x*l.y);
  //sigma_na2.x /= (l.x*l.y);
  //sigma_na2.y /= (l.x*l.y);
  
  //cout << confid_int << " " << sigma.x << " " << sigma.y << " " 
  //<< sigma_nb.x << " " << sigma_nb.y << " " << sigma_na2.x << " " << sigma_na2.y << endl;
  
  cout << confid_int << " " << X << " " << sigma.x << " " << sigma.y << endl;
  
  return 0;
}
