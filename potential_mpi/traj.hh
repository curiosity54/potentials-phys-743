#ifndef TRAJ_HH
#define TRAJ_HH

/* -------------------------------------------------------------------------- */
#include<vector>
#include<string>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include <limits>
#include<cmath>
#include <sstream>
#include<mpi.h>
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
};

class Kvectors : public traj {

int nk, n_side, nmax, loc_nmax;
double sg, kcut, invcell, kcut2;
//std::vector<double> kmod;
//std::vector<coordinates> kvec;
std::vector<double> potential;
string filename;
std::stringstream outfile;
write_file writeobj;
std::vector<double> kcoeffs; 
MPI_Comm communicator;
int rank, size, loc_nk;
std::vector<double> loc_kmod;
std::vector<coordinates> loc_kvec; 
int nat;
std::vector<double> loc_potential; 

public:

Kvectors(string infile, int N, double sigma=0.3, int nkmax=100000, MPI_Comm comm=MPI_COMM_NULL) : traj(infile), communicator(comm) {
sg = sigma;
kcut = 2.0*M_PI/sg;
kcut2= kcut*kcut;
n_side=N; //maximum number of k-vectors expected along the direction of a given cell vector
nmax = nkmax; //maximum number of k-vectors expected
invcell = get_invcell_dim();
nk =0; loc_nk=0;
outfile << "Potential_" << infile << ".npy";
writeobj.set_outfile_name(outfile.str());

MPI_Comm_size(communicator, &size);
MPI_Comm_rank(communicator, &rank); 

//if (nmax%size==0) loc_nmax= nmax/size;
//else loc_nmax = nmax/size; 

nat =get_nat();
potential.resize(nat,0);
//cout<<"Potential size in constructor  =  " <<potential.size(); 
loc_potential.resize(nat,0);
//cout<<"from constructor locnmax    "<<loc_nmax;
} ;

void gen_kvec();
double prefunc_k (double k);
void fourier_to_real();
void precalc_kcoeffs();
} ;

#endif
