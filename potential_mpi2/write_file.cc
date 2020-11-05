#include "write_file.hh"
#include<string>
#include<fstream>
#include <sstream>
#include<vector>
#include<iostream>
#include<mpi.h>


write_file::write_file(string nam, MPI_Comm comm) : communicator(comm) {filename=nam;} //cout<<"write obj declared"<<filename; }

void write_file::set_communicator(MPI_Comm comm){communicator= comm;} 
void write_file::set_outfile_name(string nam) {filename = nam;}

void write_file::writep(const vector<double> &vect, MPI_Offset offset) {
//std::ofstream fout;
//fout.open(filename);
//for(size_t i=0; i<vect.size(); i++)
//{fout<<vect[i]<<std::endl;
//}
//
MPI_File fh;
MPI_Status status;
MPI_File_open(communicator, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
int nat=vect.size();
MPI_Datatype filetype;
MPI_Comm_rank(communicator, &rank);
MPI_Comm_size(communicator, &size);
int write_size = 10;
MPI_Offset disp;
disp = rank*sizeof(double)*write_size;
disp += offset*sizeof(double)*write_size*nat;
//cout<<"    " <<rank<< "      "<<disp<<"\n";
MPI_Type_vector(nat/write_size, write_size, size*write_size, MPI_DOUBLE, &filetype);
MPI_Type_commit(&filetype);

disp = sizeof(double);
MPI_File_set_view(fh, disp, MPI_DOUBLE, filetype,"native", MPI_INFO_NULL);
MPI_File_write_all(fh, vect.data(), nat , MPI_DOUBLE, &status);

MPI_File_close(&fh);
//cout<<"\n Written successfully";
}

