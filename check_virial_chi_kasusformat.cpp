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

#define READ_NAREF

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

/***************************/

int main(int argc, char** argv) {
    string path = argv[1];
    
    /************ PARSE JSON ****************/
    string jsonstr,line;
    
    // Read JSON
    ifstream f;
    string jsonpath = path+"/in.json";
    f.open(jsonpath.c_str());
    
    while (getline(f, line)) {
        //if(line.find("## INPUT JSON ##") != string::npos)
        //  continue;
        //if(line.find("## END OF INPUT JSON ##") != string::npos)
        //  break;
        
        jsonstr += line;
        jsonstr += '\n';
    }
    
    f.close();
    cout << "Parsing JSON block..." << endl;
    //cout <<jsonstr << endl;
    
    //remove #
    replaceAll(jsonstr, "#", "");
    
    Document d;
    d.Parse<kParseCommentsFlag>(jsonstr.c_str());
    
    // parse require
    int Nx = 32; //d["N"][0].GetInt();
    int Ny = 32; //d["N"][1].GetInt();
    int N = Nx*Ny;
    double rc = 2.5; //d["nb_rc"].GetDouble();
    //double eps = d["epsilon"].GetDouble()[0];
    double hchi = d["hchi"].GetDouble();
    Vec2d l,lh,a;
    
    cout << rc << " " << hchi << endl;
    
    //*********************** Read reference ***********************//
    
    vector<Particle> pref;
    double ignore;
    double lxl,lxr,lyl,lyr;
    int pcount = 0;
    
    f.open(path+"/naref.txt");
    
    getline(f, line); // # N X ibin lx_l lx_r ly_l ly_r mid post pre
    
    getline(f, line); // #
    replaceAll(line, "# ", "");
    
    //stringstream iss(line);
    //iss >> N >> ignore >> ignore >> lxl >> lxr >> lyl >> lyr >> ignore >> ignore >> ignore;
    
    vector<string> items = split(line, ' ');
    N = stoi(items[0]);
    l.x = stod(items[4])-stod(items[3]);
    l.y = stod(items[6])-stod(items[5]);
    lh = l*0.5;
    
    getline(f, line); // # pid type x y chi
    
    cout << "ref -> " << N << " " << " " << l.x << " " << l.y << endl;
    pref.resize(N);
    while (getline(f, line)) {
        stringstream iss(line);
        iss >> pref[pcount].pid
        >> ignore
        >> pref[pcount].r.x
        >> pref[pcount].r.y
        >> pref[pcount].chi;
        
        pcount++;
    }
    
    f.close();
    
    /********************** Init NonAffine ******************/
    NonAffine nf;
    nf.init(Nx,Ny,l,hchi);
    nf.set_ref(pref);
    
    vector<Particle> particles;
    particles.resize(N);
    
    NonBondedForces nb;
    
    vector<Particle> particles_acc;
    particles_acc.resize(N);
    
    /*********************** Read Particles ******************/
    string header_str;
    int en_id = 10;
    
    int confid_begin = stoi(argv[2]);
    int confid_end   = stoi(argv[3]);
    int nens         = stoi(argv[4]);
    
    //
    cout << "--------- Lets begin ----------" << endl;
    char fname[128];
    string confpath;
    
    for(int confid=confid_begin; confid<=confid_end; confid++) {
        
        // averaging stuffs
        for(auto &ip:particles_acc) {
            ip.pid = 0;
            ip.r = Vec2d(0,0);
            ip.chi = 0;
            ip.pe = 0;
            for(int j=0; j<4; j++) ip.virial_tensor[j] = 0.0;
            for(int j=0; j<4; j++) ip.virial_na2[j] = 0.0;
        }

        
        
        for(int ensid=0; ensid<=nens; ensid++) {
            snprintf(fname,128,"/cbin_%04d_%05d.conf",confid,ensid);
            
            confpath = path+string(fname);
            f.open(confpath, ios::in);
            
            if(!f.is_open() ) {
                cerr << "# cant open file " << confpath << endl;
                break;
            }
            
            
            int bin_id;
            double X_bin, bin_id_double;
            
            //cout << confpath << endl;
            
            getline(f, line);
            
            double tmp_d;
            int ign;
            
            getline(f, line);
            replaceAll(line, "#", "");
            
            stringstream iss1(line);
            iss1 >> tmp_d >> tmp_d >> bin_id_double >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore;
            bin_id = int(bin_id_double);
            
            getline(f, line);
            
            //cout << "header -> " << bin_id << " " << X_bin << " " << l.x << " " << l.y << endl;
            
            pcount = 0;
            
            while (getline(f, line)) {
                stringstream iss(line);
                iss >> particles[pcount].pid
                >> ign
                >> particles[pcount].r.x
                >> particles[pcount].r.y
                >> particles[pcount].chi;
                
                pcount++;
                
            }
            
            //for(auto &i:particles)
            //    cout << i.pid << " " << i.r.x << " " << i.r.y << " " << endl;
            
            f.close();
            
            double X = nf.calc_chi(particles);
            
            nb.init(rc, l, N, particles);
            nb.calc();
            nb.get_results(particles);
            
            ofstream of(confpath+".wstress", ios::out);
            of << header_str << endl;
            for(auto &ip:particles) {
                of << ip.pid << " "
                << ip.r.x << " "
                << ip.r.y << " "
                << ip.chi << "\t"
                << ip.virial_tensor[0] + ip.virial_na2[0] << " "
                << ip.virial_tensor[3] + ip.virial_na2[1] << " "
                << ip.virial_na2[0] << " "
                << ip.virial_na2[1] << endl;
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
            
            // Add to particles_acc
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
                particles_acc[i].virial_tensor[1] += particles[i].virial_tensor[3];
                particles_acc[i].virial_na2[0] += particles[i].virial_na2[0];
                particles_acc[i].virial_na2[1] += particles[i].virial_na2[1];
            }

            
            cout << confid << " " << ensid << " " << X << " " << sigma.x << " " << sigma.y << endl;
        }
        
        // Save average data
        ofstream ofav;
        
        char str_bin_id[64];
        snprintf(str_bin_id,64,"/av_%04d.txt",confid);
        
        ofav.open(path+str_bin_id);
        int av_over_bins_counter = nens + 1;
        
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
    
    return 0;
}
