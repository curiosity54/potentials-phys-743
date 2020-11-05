/* -------------------------------------------------------------------------- */
#include "traj.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <iostream>
#include <string> 
#include <stdio.h>

using namespace std; 
using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

#define N 50
#define NMAX 100000
#define SIGMA 0.3
int main()
{ 

string filename="traj_20.xyz";
//traj t1("traj.xyz");
//t1.print_traj();
//t1.center_coords();
//cout<<t1.get_invcell_dim(); 
Kvectors k1(filename, N, SIGMA, NMAX);

auto t1 = clk::now();
k1.gen_kvec();
k1.fourier_to_real();
second elapsed = clk::now() - t1;


std::printf("\n wall clock time (chrono)        = %.4gs\n", elapsed.count());

return 0;

}

