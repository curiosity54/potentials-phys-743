#ifndef WRITE_FILE_HH
#define WRITE_FILE_HH

#include<string>
#include<fstream>
#include <sstream>
#include<vector>
using namespace std;
class write_file {
string filename;
public :
write_file(string nam="Potential.npy");
void writep(const vector<double> &vect);
};


#endif
