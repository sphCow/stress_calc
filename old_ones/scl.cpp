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

int main(int argc, char** argv) {

	int start_from_bin = -1;
	
	string path = argv[1];
	bool export_per_particle_data = false;
	
	//if(path[pat])
	
	cout << "path --> " << path << endl;
	
	
	
	string line;
	string jsonstr;
	
	
	cout << "Reading JSON block..." << endl;
	// parse json and get parameters
	ifstream f;
	f.open(path.c_str());
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
	cout << "Parsing JSON block..." << endl;
	
	//remove #
	replaceAll(jsonstr, "#", "");
	
	cout << jsonstr << endl;
	
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
	
	cout << N << " " << rc << " " << hchi << endl;
	
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
		
	string refpath = path+"ref";
	ofstream fref(refpath);
	for(int i=0; i<N; i++) {
		fref << pref[i].pid << " " << pref[i].r.x << " " << pref[i].r.y << endl;
	}
	
	//*******************END of Generate reference ***********************//
	
	
	//******************* Init NonAffine ********************* //
	NonAffine nf;
    nf.init(Nx,Ny,l,hchi);
    nf.set_ref(pref);
	
	// Parse each block and compute stress 
	cout << "Parsing data blocks..." << endl;
	f.open(path.c_str());
	lno = 0;
	int pcount = 0;
	size_t target_lno = 1e10;
	int bin_id = -1;
	int bin_id_prev = -1;
	double ignore,X_bin;
	string header_str;
	int av_over_bins_counter;
	
	vector<Particle> particles;
	vector<Particle> particles_acc;
	
	particles.resize(N);
	particles_acc.resize(N);
	cout << "Particles array resized to " << N << endl;
	
	NonBondedForces nb;
	ofstream of, of_thermo;
	
	//of.open(path+".wstress");
	of_thermo.open(path+".wstress.thermo");
	of_thermo << "#bin_id  X_bin <X> eff_stress stress_xx stress_yy stress_xy" << endl;
	
	while (getline(f, line)) {
		size_t found1 = line.find("# NEW_BIN");
		
		if (found1 != string::npos) {
			target_lno = lno;
			pcount = 0;
			
			header_str = line;
			replaceAll(line, "INTERMIDIATE", "");
			replaceAll(line, "#", "");
			replaceAll(line, "NEW_BIN", "");
			double bin_id_double;
			
			stringstream iss(line);
			iss >> bin_id_double >> ignore >> X_bin >> ignore >> l.x >> l.y;
			bin_id = bin_id_double*2;
			
			// Handle average over bin stuffs
			if(bin_id>=start_from_bin &&  bin_id_prev != bin_id) {
				cout << "NEW_BIN_ID " << bin_id << endl;
				
				//********* Export quantities averaged over bins ***************//
				if(bin_id>0) {
					stringstream ss;
					replaceAll(path,".bconfig","");
					ss << path << "_" << bin_id_prev << ".avbin";
					cout << "File opened --> " << ss.str() << endl;
					
					ofstream f(ss.str(), ios::out) ;
					double c = static_cast<double>(av_over_bins_counter);
					f << header_str << " c=" << c << endl;
					
					for(auto &ip: particles_acc) {
						f << ip.pid/c << " " 
							<< ip.r.x/c << " " 
							<< ip.r.y/c << " " 
							<< ip.chi/c << " " 
							<< ip.virial_tensor[0]/c << " " 
							<< ip.virial_tensor[3]/c << " "
							<< ip.virial_tensor[1]/c << " "
							<< ip.pe/c << endl;
					}
					
					f.close();
				}
				
				//Reset accumualtors
				av_over_bins_counter = 0;
				for(auto &ip:particles_acc) {
					ip.pid = 0;
					ip.r = Vec2d(0,0);
					ip.chi = 0;
					ip.pe = 0;
					for(int j=0; j<4; j++) ip.virial_tensor[j] = 0.0;
				}
				
				//Set new bin id
				bin_id_prev = bin_id;
				
			}
			
			av_over_bins_counter++;
			cout << "header --> " << bin_id << " " << av_over_bins_counter << "\r";
			cout.flush();
		}
		
		if(bin_id >= start_from_bin) {
			if(lno>target_lno && pcount<N) {
				//cout << lno << " " << pcount << " --> " << line << endl;
				
				stringstream iss(line);
				iss >> particles[pcount].pid 
						>> particles[pcount].r.x
						>> particles[pcount].r.y
						>> particles[pcount].chi;
				
				pcount++;
				
			}

			//done reading particles in this bin
			else if(pcount==N) {
				pcount = 0;
				
				//Do whatever you want to do with the particles array here
				nb.init(rc, l, N, particles);
				nb.calc();
//				nb.get_results(particles);
                
                /************* NonAffine stuffs to be added here **************/
                nf.calc_chi(particles);
                
                for(auto &ip:particles) {
                  if(fabs(ip.chi-ip.chi2)>0.0001) cerr << "! the guys are not innocent !" << endl;
                  ip.virial_tensor[0] = ip.r.x*ip.f.x;
                  ip.virial_tensor[1] = ip.r.x*ip.f.y;
                  ip.virial_tensor[2] = ip.r.y*ip.f.x;
                  ip.virial_tensor[3] = ip.r.y*ip.f.y;
                }
				
				//********* Export Per particle data ***************//
				if (export_per_particle_data){
					of.open(path + "_" + to_string(bin_id) + ".wstress");
					of << header_str << endl;
					for(auto &i:particles){
						of << i.pid << " " 
							   << i.r.x << " "
							   << i.r.y << " "
							   << i.chi << " "
                               << i.virial_tensor[0] << " " 
                               << i.virial_tensor[3] << " "
                               << i.virial_tensor[1] << endl;
								 
					}
					of.close();
				}
				
				//****************** Global thermo *****************//
				double X = 0.0;
				//double stress = 0.0;
				double stress_xx = 0.0;
				double stress_yy = 0.0;
				double stress_xy = 0.0;
				
				/*
				i .virial_tensor[0] += ff* rr.x*rr.x; //xx   *
				i.virial_tensor[1] += ff* rr.y*rr.x; //xy
				i.virial_tensor[2] += ff* rr.x*rr.y; //yx
				i.virial_tensor[3] += ff* rr.y*rr.y; //yy
				 */
				
				for(auto &i:particles) {
					X += i.chi;
					//stress += (i.virial_tensor[3] - i.virial_tensor[0]);
					stress_xx += i.virial_tensor[0];
					stress_xy += i.virial_tensor[1];
					stress_yy += i.virial_tensor[3];
				}
				
				X /= double(N);
				stress_xx /= double(l.x*l.y);
				stress_xy /= double(l.x*l.y);
				stress_yy /= double(l.x*l.y);
				double eff_stress = 0.5*(stress_xx-stress_yy);
				
				of_thermo << bin_id << " " << X_bin << " " << X << " " 
				<< eff_stress << " " << stress_xx << " " << stress_yy << " " << stress_xy << endl;
				//************** End Global thermo *****************//
				
				//Accumulate
				for(int i=0; i<N; i++) {
					
					//PBC 
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
					
					for(int j=0; j<4; j++) 
						particles_acc[i].virial_tensor[j] += particles[i].virial_tensor[j];
					
					particles_acc[i].pid += particles[i].pid;
					particles_acc[i].pe += particles[i].pe;
				}
				
			}
	} else {
		cout << "bin_id =" << bin_id << " --> SKIPPED" << endl;
	}
		
		lno++;
	}
	
	
	of.close();
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return 0;
}
