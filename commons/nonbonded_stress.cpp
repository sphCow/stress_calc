#include "nonbonded_stress.h"

void NonBondedForces::init(double rc_, Vec2d l_, int N_, vector<Particle>& particles_) {
	rc = rc_;
	l = l_;
	N = N_;
	particles = particles_;
	
	lh = l*0.5;
	rc2 = rc*rc;
	
  ncell_x = floor(l.x/rc);  //2;
  ncell_y = floor(l.y/rc);  //2;
  ncell = ncell_x * ncell_y;

  lcell.x = l.x / double(ncell_x);
  lcell.y = l.y / double(ncell_y);

  //cout << "cells " << ncell_x << "x" << ncell_y << endl;

  //cell resize
  clen.resize(ncell);
  cbegin.resize(ncell);

  c0 = 0.0;
  c2 = 0.0;
  c6 = 0.0;
}

void NonBondedForces::update_cell_list() {
    for(int i=0; i<ncell; i++) {
        clen[i] = 0;
    }

    for(int i=0; i<N; i++) {
        int cx = ((particles[i].r.x + lh.x) / lcell.x);
        int cy = ((particles[i].r.y + lh.y) / lcell.y);
        int cid = (cx) + (cy) * (ncell_x);
        //cout << cx << " " << cy << " " << cid << endl;
        particles[i].cid = cid;
        clen[cid]++;
    }

    sort(particles.begin(), particles.end(), sort_cid);

    int sum=0;
    for(int i=0; i<ncell; i++) {
        cbegin[i] = sum;
        sum += clen[i];
    }
}

void NonBondedForces::get_pair_LJ(Particle &i, Particle &j) {
    Vec2d rr = i.r - j.r;
    pbc(rr);

    double r2 = rr.x*rr.x + rr.y*rr.y;

    if(r2<rc2) {
        double rri = 1.0/r2; //r^-2

        double r14i = pow(rri, 7.0);
        double r8i = pow(rri, 4.0);
        double r12i = pow(rri, 6.0);
        double r6i = pow(rri, 3.0);
        double ff = 48*r14i - 24*r8i - 8*c2 - 24*r2*r2*c6;

        i.f.x += 0.5*rr.x*ff;
        i.f.y += 0.5*rr.y*ff;
        //j.f -= rr*ff;
				
		double this_pe = 4.0*(r12i - r6i);
				
        pe += this_pe;
		i.pe += this_pe;

        //virial_trace += ff*r2;
        i.virial_tensor[0] += ff* rr.x*rr.x; //xx
        i.virial_tensor[1] += ff* rr.y*rr.x; //xy
        i.virial_tensor[2] += ff* rr.x*rr.y; //yx
        i.virial_tensor[3] += ff* rr.y*rr.y; //yy

        //double phi = 32.0*rri*(21.0*r14i - 6.0*r8i);
        //double phi = 96.0*rri*(7.0*r14i - 2.0*r8i + r2*r2*c6);
        //double phi = rri*(672.0*r14i - 192*r8i + 96*c6*r2*r2);
        //born_xx += xr*xr*xr*xr*phi; // 1111
        //born_xy += xr*xr*yr*yr*phi; // 1122
        //born_yx += xr*xr*yr*yr*phi; // 1122
        //born_yy += yr*yr*yr*yr*phi;
    }
}

void NonBondedForces::calc() {
    update_cell_list();
		
		for(auto &ip:particles) {
			ip.virial_tensor[0] = 0.0;
			ip.virial_tensor[1] = 0.0;
			ip.virial_tensor[2] = 0.0;
			ip.virial_tensor[3] = 0.0;
            ip.f.x=0.0;
            ip.f.y=0.0;
		}
		
    int i,j;

    //#pragma omp parallel for shared(en,virial_sum,cells)
    for (int cell_id=0; cell_id<ncell; cell_id++) { //loop over all cells ic in block b

        int iy = int(cell_id/ncell_x);
        int ix = cell_id%ncell_x;
				
				int dir[8];
				
        dir[0] = ((ix + 1) % ncell_x + ((iy + ncell_y) % ncell_y) * ncell_x); //EAST
				dir[1] = ((ix + 1) % ncell_x + ((iy + ncell_y + 1) % ncell_y) * ncell_x); //NE
				dir[2] = ((ix + ncell_x) % ncell_x + ((iy + ncell_y + 1) % ncell_y) * ncell_x); //N
				dir[3] = ((ix + ncell_x - 1) % ncell_x + ((iy + ncell_y + 1) % ncell_y) * ncell_x); //NW
				dir[4] = ((ix + ncell_x - 1) % ncell_x + ((iy + ncell_y) % ncell_y) * ncell_x); //W
				dir[5] = ((ix + ncell_x - 1) % ncell_x + ((iy + ncell_y - 1) % ncell_y) * ncell_x); //SW
				dir[6] = ((ix + ncell_x) % ncell_x + ((iy + ncell_y - 1) % ncell_y) * ncell_x); //S
        dir[7] = ((ix + 1) % ncell_x + ((iy + ncell_y - 1) % ncell_y) * ncell_x); //SE
        //(ix + iy * ncell_x); //base_cell //own cell */

				// Neighbour cell
        for (i=cbegin[cell_id]; i<cbegin[cell_id]+clen[cell_id]; i++ ) {       // loop over particles in i cell
          for(int d=0;d<8;d++) {
						for (j = cbegin[dir[d]]; j<cbegin[dir[d]]+clen[dir[d]]; j++ ) {
	                get_pair_LJ(particles[i], particles[j]);
	            }
					}
        }

        // Own cell
        for(int i = cbegin[cell_id]; i<cbegin[cell_id]+clen[cell_id]; i++ ) {
					for(int j = cbegin[cell_id]; j<cbegin[cell_id]+clen[cell_id]; j++) {
                if(i != j) get_pair_LJ(particles[i],particles[j]);
            }
        }
    
			
		}
		
}
