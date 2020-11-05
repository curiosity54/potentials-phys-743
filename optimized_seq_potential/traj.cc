#include "traj.hh"

#include<vector>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<sstream>
using namespace std; 


traj::traj(string file):filename(file)
{
  fstream f;
  coordinates cc;
 f.open(filename, fstream::in);
  f>>this->natoms;
  cout<<"\n nat "<<natoms;
  f>>this->cell;
 coords.resize(natoms);

 // cout<<"\n cell "<<cell;
  for (int i=0; i<natoms; i++)
{
 f>>cc.x>>cc.y>>cc.z;
coords[i]=cc;
}
//cout<<"Num of coords read "<<coords.size();
this->volume= std::pow(cell,3);
centered_coords.resize(natoms*(natoms-1));
}


void traj::print_traj() {
  for (int i=0; i<natoms; i++)
   {cout<<"\n"<<this->coords[i].x <<this->coords[i].y<<this->coords[i].z;
       }
}


double traj::get_invcell_dim() {
return 2*M_PI/cell;
}



std::vector<coordinates> traj::center_coords() {
coordinates diff(0,0,0);
int idiff=0;
for (int i=0; i<natoms; i++) {
for (int j=0; j<natoms; j++) {
if (i!=j)
{diff = this->coords[j] - this->coords[i];
centered_coords[idiff]=diff; idiff++; } 
}}
//  DEBUG  cout<<"\n"<<centered_coords.size()<<"difference coords created \n";
//return centered_coords;
//for(int i =0; i<centered_coords.size(); i++) {        
//cout<<"\n"<<i << centered_coords[i].x << centered_coords[i].y <<centered_coords[i].z;
//}
return centered_coords;
}

void Kvectors::gen_kvec() {
int ksize=600000;
kmod.resize(ksize); 
kvec.resize(ksize);
    //compute the k-vectors on a semi-sphere of radius kcut (origin excluded) oriented along b1
    //the half-line i1=0,i2=0,i3>0
    for(int i3=1; i3<=n_side; i3++){
      //compute the k-vector
      double kx = 0.;
      double ky = 0.;
      double kz = i3*invcell;
      double k2 = kx*kx + ky*ky + kz*kz;
      if (k2<=kcut2){
         kmod[nk]=std::sqrt(k2);
         kvec[nk]=coordinates(kx,ky,kz);
         nk += 1;
      }
    }

    //the half-plane i1=0,i2>0
   for(int i2=1; i2<=n_side; i2++){
      for(int i3=0; i3<=2*n_side; i3++){
        int ivec3 = i3 - n_side;
        //compute the k-vector
        double kx = 0.;
        double ky = i2*invcell;
        double kz = ivec3*invcell;
        double k2 = kx*kx + ky*ky + kz*kz;
        if (k2<=kcut2){
         // kmod.push_back(std::sqrt(k2));kvec.push_back(coordinates(kx,ky,kz));
         kmod[nk]=std::sqrt(k2);
         kvec[nk]=coordinates(kx,ky,kz);
          nk += 1;
        }
      }
    }
    //the half-space i1>0
  for(int i1=1; i1<=n_side; i1++){
      for(int i2=0; i2<=2*n_side; i2++){
        int ivec2 = i2 - n_side;
        for(int i3=0; i3<=2*n_side; i3++){
          int ivec3 = i3 - n_side;
          //compute the k-vector
          double kx = i1*invcell;
          double ky = ivec2*invcell;
          double kz = ivec3*invcell;
          double k2 = kx*kx + ky*ky + kz*kz;
          if (k2<=kcut2){
          kmod[nk]=std::sqrt(k2);
          kvec[nk]=coordinates(kx,ky,kz);
          //kmod.push_back(std::sqrt(k2));
          //kvec.push_back(coordinates(kx,ky,kz));
          nk += 1;}
             }
           }
         }
printf("\n %d number of kvecs generated", nk);
//cout<<"\n" <<"dim of moduli"<< kmod.size();
//cout<<"\n" <<"dim of vec"<< kvec.size();
}
double Kvectors::prefunc_k (double k) {
return std::exp(-0.5*(sg*k*sg*k)) / (k*k) * std::sin(k);
}

void Kvectors::precalc_kcoeffs() {
kcoeffs.resize(nk);
for(int i=0;i<nk; i++) 
{kcoeffs[i]= prefunc_k(kmod[i]); 
}
}

void Kvectors::fourier_to_real()
{double real, arg, kx, ky, kz, coordx, coordy, coordz, sum=0.;
std::vector<coordinates> centered_coords = center_coords();
coordinates diff;
//int nat =get_nat();
//potential.resize(nat,0); 
precalc_kcoeffs();
for (int i=0; i<nat; i++) {
sum=0.;
    for (int iG=0; iG<nk; iG++) {
           kx = kvec[iG].x;
           ky = kvec[iG].y;
           kz = kvec[iG].z;
            real=0.;
       for(int j=0; j<nat-1; j++) {
            coordinates diff=  centered_coords[i*nat+j-i];
            coordx= diff.x;
            coordy= diff.y;
            coordz= diff.z;
            arg = kx*coordx + ky*coordy + kz*coordz;
            real+= std::cos(arg);
         }
      sum+= real * kcoeffs[iG];
     }
   sum = sum /get_vol(); // * 16* pow(M_PI,2) 
   potential[i]=sum;
   cout<<"\n pot"<<i<<sum;}
   
cout<<"\n Potential real space dime"<<potential.size();
writeobj.writep(potential ) ;
}
