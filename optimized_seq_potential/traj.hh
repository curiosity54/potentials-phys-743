#ifndef TRAJ_HH
#define TRAJ_HH

/* -------------------------------------------------------------------------- */
#include<vector>
#include<string>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include <sstream>
/* -------------------------------------------------------------------------- */
#include "coordinates.hh"
#include "write_file.hh"
/* -------------------------------------------------------------------------- */
using namespace std; 
class traj {

string filename;
double cell, volume;
int natoms;

std::vector<coordinates> coords;
std::vector<coordinates> centered_coords;

public:

traj(string file);
void print_traj();
double get_invcell_dim();

double get_vol() const {return volume;}

std::vector<coordinates> center_coords();

int get_nat() const {return natoms;}

} ;


class Kvectors : public traj {

int nk, n_side, nmax, nat;
double sg, kcut, invcell, kcut2;
std::vector<double> kmod;
std::vector<coordinates> kvec;
std::vector<double> potential;
string filename;
std::stringstream outfile;
write_file writeobj;
std::vector<double> kcoeffs; 
public:

Kvectors(string infile, int N, double sigma=0.3, int nkmax=100000) : traj(infile)  {
sg = sigma;
kcut = 2.0*M_PI/sg;
kcut2= kcut*kcut;
n_side=N; //maximum number of k-vectors expected along the direction of a given cell vector
nmax = nkmax; //maximum number of k-vectors expected
invcell = get_invcell_dim();
nk =0;
outfile << "Potential_" << infile << ".npy";
writeobj.set_outfile_name(outfile.str());
nat = get_nat(); 
potential.resize(nat,0); 
} ;
//
void gen_kvec();
double prefunc_k (double k);
void fourier_to_real();
void precalc_kcoeffs();
} ;
//
//

#endif
