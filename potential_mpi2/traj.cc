#include "traj.hh"

#include<vector>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include<string>
#include<cmath>
#include<sstream>
#include<mpi.h>
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

//  Read a chunk of atoms 
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



std::vector<coordinates> traj::center_coords(int begin, int end) {
coordinates diff(0,0,0);
int idiff=0;
//cout<<"\n begin, end from center_coords  " <<begin<<end;  
for (int i=begin; i<end; i++) {
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
loc_kmod.resize(ksize); 
loc_kvec.resize(ksize);
int start_i1, start_i2, start_i3, end_i1, end_i2, end_i3; 
start_i3 =  n_side/(new_size)*new_rank+1; 
end_i3 = (new_rank<new_size-1)? n_side/new_size*(new_rank+1) : n_side;     
//cout<< "rank ,start_i3   "<< rank<<start_i3; 
//cout<< "rank ,end_i3   "<< rank<<end_i3;
//compute the k-vectors on a semi-sphere of radius kcut (origin excluded) oriented along b1
    //the half-line i1=0,i2=0,i3>0
    for(int i3=start_i3; i3<=end_i3; i3++){
      //compute the k-vector
      double kx = 0.;
      double ky = 0.;
      double kz = i3*invcell;
      double k2 = kx*kx + ky*ky + kz*kz;
      if (k2<=kcut2){
         loc_kmod[loc_nk] = std::sqrt(k2);//loc_kmod.push_back(std::sqrt(k2));
         loc_kvec[loc_nk] = coordinates(kx,ky,kz); //loc_kvec.push_back(coordinates(kx,ky,kz));
         loc_nk += 1;
      }
    }

start_i2 = n_side/(new_size)*new_rank+1; //n_side/size * rank+1;
end_i2 = (new_rank<new_size-1)? n_side/new_size*(new_rank+1) : n_side; //(rank<size-1)? n_side/size*(rank+1) : n_side;
    //the half-plane i1=0,i2>0
   for(int i2=start_i2; i2<=end_i2; i2++){
      for(int i3=0; i3<=2*n_side; i3++){
        int ivec3 = i3 - n_side;
        //compute the k-vector
        double kx = 0.;
        double ky = i2*invcell;
        double kz = ivec3*invcell;
        double k2 = kx*kx + ky*ky + kz*kz;
        if (k2<=kcut2){
          loc_kmod[loc_nk] = std::sqrt(k2);//loc_kmod.push_back(std::sqrt(k2));
          loc_kvec[loc_nk] = coordinates(kx,ky,kz);//loc_kvec.push_back(coordinates(kx,ky,kz));
          loc_nk += 1;
        }
      }
    }

start_i1 = n_side/(new_size)*new_rank+1;//n_side/size * rank+1;
end_i1= (new_rank<new_size-1)? n_side/new_size*(new_rank+1) : n_side;//(rank<size-1)? n_side/size*(rank+1) : n_side;
    //the half-space i1>0
  for(int i1=start_i1; i1<=end_i1; i1++){
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

          loc_kmod[loc_nk] = std::sqrt(k2);//loc_kmod.push_back(std::sqrt(k2));
          loc_kvec[loc_nk] = coordinates(kx,ky,kz);//loc_kvec.push_back(coordinates(kx,ky,kz));
          loc_nk += 1;}
             }
           }
         }

MPI_Allreduce(&loc_nk, &nk, 1, MPI_INT, MPI_SUM, newcomm); 

//if (new_rank ==0) {printf("\n %d  total number of kvecs generated",  nk);}
//{printf("\n %d  %d   %d world rank, rank, total number of kvecs generated", rank, new_rank,  nk);}
}

double Kvectors::prefunc_k (double k) {
return std::exp(-0.5*(sg*k*sg*k)) / (k*k) * std::sin(k);
}

void Kvectors::precalc_kcoeffs() {
kcoeffs.resize(loc_nk);
for(int i=0;i<loc_nk; i++) 
{kcoeffs[i]= prefunc_k(loc_kmod[i]); 
}
//cout<<"\n kcoeff shape   "<< kcoeffs.size(); 
}

void Kvectors::fourier_to_real()
{double real, arg, kx, ky, kz, coordx, coordy, coordz, sum=0.;
std::vector<coordinates> centered_coords = center_coords(begin, end);
coordinates diff;
int idx=0;
precalc_kcoeffs(); int idiff=0;
for (int i=begin; i<end; i++) {
sum=0.;
    for (int iG=0; iG<loc_nk; iG++) { idiff=idx*(nat-1); 
           kx = loc_kvec[iG].x;
           ky = loc_kvec[iG].y;
           kz = loc_kvec[iG].z;
           real=0.;
       for(int j=0; j<nat-1; j++) {
            
            coordinates diff=  centered_coords[idiff];
            coordx= diff.x;
            coordy= diff.y;
            coordz= diff.z;
            arg = kx*coordx + ky*coordy + kz*coordz;
            real+= std::cos(arg);idiff++;
         }
      sum+= real * kcoeffs[iG];
     }
   sum = sum /get_vol(); // * 16* pow(M_PI,2) 
   loc_potential[idx]=sum;
   idx++;
//cout<<"\n pot from rank "<< rank <<"  " <<i<<sum;
   }
MPI_Allreduce(loc_potential.data(), potential.data(), len, MPI_DOUBLE, MPI_SUM, newcomm);   
//cout<<"\n \n Potential real space dime"<<potential.size();
if (new_rank ==0)
{//cout<<"idx from rank   "<<idx<<new_rank; //cout<<"\n \n Potential real space dime"<<potential.size();
for(int i=0; i<idx; i++)
cout<<"\n pot"<<i+begin<<"   " <<potential[i]; 
}
MPI_Offset offs; offs = begin; 
writeobj.set_communicator(newcomm);
writeobj.writep(potential,offs) ;
MPI_Comm_free(&newcomm); 
}
