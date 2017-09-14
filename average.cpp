#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
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

std::vector<std::string> split(const std::string &text, char sep) {
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        if (end != start) {
            tokens.push_back(text.substr(start, end - start));
        }
        start = end + 1;
    }
    if (end != start) {
        tokens.push_back(text.substr(start));
    }
    return tokens;
}

/********************* read particles *********************/
int read_particles(string path, vector<Particle> &particles, long &bin_id) {
  ifstream f;
  f.open(path,ios::in);
  
  if(!f.good()) {
    cerr << "cant find file " << path << endl;
  }
  
  int pcount=0;
  double X_bin, bin_id_double;
  string header_str,ignore;
  Vec2d l;
  
  string line;
  getline(f, line);
  header_str = line;
  
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
        >> particles[pcount].chi
        >> particles[pcount].virial_tensor[0]
        >> particles[pcount].virial_tensor[1]
        >> particles[pcount].virial_na2[0]
        >> particles[pcount].virial_na2[1];
    
    pcount++;
    
  }
  
  f.close();
}

int main(int argc, char** argv) {
  
  string path = argv[1];
  
  string line;
  string jsonstr;
  
  //************************* Parse JSON block **************************//
  cout << "Reading JSON block..." << endl;
  // parse json and get parameters
  ifstream f;
  f.open(path+"/c000000");
  bool json_block = false;
  size_t lno = 0;
  
  while (getline(f, line)) {
    size_t found1 = line.find("## END OF INPUT JSON ##");
    if (found1 != string::npos) {
      json_block = false;
      break;
    }
    
    if(json_block) {
      jsonstr += line;
      jsonstr += '\n';
    }
    
    size_t found = line.find("## INPUT JSON ##");
    
    if (found != string::npos) {
      json_block = true;
    }
    
    lno++;
  }
  
  f.close();
  
  replaceAll(jsonstr, "#", "");
  Document d;
  d.Parse<kParseCommentsFlag>(jsonstr.c_str());
  
  // parse required fields
  int Nx = d["N"][0].GetInt();
  int Ny = d["N"][1].GetInt();
  int N = Nx*Ny;
  double rc = d["nb_rc"].GetDouble();
  double rho = d["rho"].GetDouble();
  double eps = d["epsilon"].GetDouble();
  double hchi = d["hchi"].GetDouble();
  Vec2d l,lh,a;
  
  cout << N << " " << rc << " " << hchi << " " << eps << endl;
  
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
  
  l.x += eps*l.x ;
  l.y -= eps*l.y ;
  
  for(int i=0; i<N; i++) {
    pref[i].r.x += eps*pref[i].r.x; // x
    pref[i].r.y -= eps*pref[i].r.y; 
  }
  
  //******************* read particle data *********************//
  vector<Particle> particles;
  vector<Particle> particles_acc;
  particles.resize(N);
  particles_acc.resize(N);
  
  long bin_id = -1;
  long prev_bin_id = -1;
  long av_over_bins_counter = 0;
  double X;
  
  for(long confid_int=1; confid_int<40800; confid_int++) {
    char confid[64];
    snprintf(confid,64,"/c%06d_wstress",confid_int);
  
    string confpath = path+string(confid);
    
    read_particles(confpath, particles, bin_id);
    
    if(bin_id != prev_bin_id) {
      cout << "New bin" << endl;
      
      if(bin_id>1) {
        // save average data
        ofstream ofav;
        
        char str_bin_id[64];
        snprintf(str_bin_id,64,"/c%06d_avbin",bin_id-1);
        
        ofav.open(path+str_bin_id);
        ofav << "# av over " << av_over_bins_counter << endl;
        
        for(auto &ip:particles_acc) {
          ofav << ip.pid/double(av_over_bins_counter) << " " 
          << ip.r.x/double(av_over_bins_counter) << " "
          << ip.r.y/double(av_over_bins_counter) << " "
          << ip.chi/double(av_over_bins_counter) << "\t"
          << ip.virial_tensor[0]/double(av_over_bins_counter) << " "
          << ip.virial_tensor[1]/double(av_over_bins_counter) << "\t" 
          << ip.virial_na2[0]/double(av_over_bins_counter) << " "
          << ip.virial_na2[1]/double(av_over_bins_counter) << endl;
        }
        
        ofav.close();
      }
      
      //reset accumulators
      av_over_bins_counter = 0;
      for(auto &ip:particles_acc) {
        ip.pid = 0;
        ip.r = Vec2d(0,0);
        ip.chi = 0;
        ip.pe = 0;
        for(int j=0; j<4; j++) ip.virial_tensor[j] = 0.0;
        for(int j=0; j<4; j++) ip.virial_na2[j] = 0.0;
      }
      
    }
    
    prev_bin_id = bin_id;
    cout << confid << " -> " << bin_id << endl;
    
    //************ Accumulate ************************//
    av_over_bins_counter += 1;
    
    for(int i=0; i<N; i++) {
      // pid 
      particles_acc[i].pid += particles[i].pid;
      particles_acc[i].chi += particles[i].chi;
      
      // r with pbc 
      if((particles[i].r.x-pref[i].r.x)>lh.x)
        particles[i].r.x = particles[i].r.x-l.x;
      
      else if((particles[i].r.x-pref[i].r.x)<-1*lh.x)
        particles[i].r.x = particles[i].r.x+l.x	;				
      
      if((particles[i].r.y-pref[i].r.y)>lh.y)
        particles[i].r.y = particles[i].r.y-l.y;
      
      else if((particles[i].r.y-pref[i].r.y)<-1*lh.y)
        particles[i].r.y = particles[i].r.y+l.y;
      
      particles_acc[i].r += particles[i].r;
      particles_acc[i].chi += particles[i].chi;
      
          
      // accumulate stress 
      particles_acc[i].virial_tensor[0] += particles[i].virial_tensor[0];
      particles_acc[i].virial_tensor[1] += particles[i].virial_tensor[1];
      particles_acc[i].virial_na2[0] += particles[i].virial_na2[0];
      particles_acc[i].virial_na2[1] += particles[i].virial_na2[1];
    }
    
    
  
  }
  
  return 0;
}
  
  
