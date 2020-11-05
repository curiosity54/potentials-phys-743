#ifndef WRITE_FILE_HH
#define WRITE_FILE_HH

#include<string>
#include<fstream>
#include <sstream>
#include<vector>
#include<mpi.h>
using namespace std;
class write_file {
string filename;
MPI_Comm communicator; 
int rank, size; 
public :
write_file(string nam="Potential.npy", MPI_Comm = MPI_COMM_WORLD);
void set_outfile_name(string nam);
void writep(const vector<double> &vect, MPI_Offset offset);
void set_communicator(MPI_Comm comm);  
};


#endif
