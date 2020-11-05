#include "write_file.hh"
#include<string>
#include<fstream>
#include <sstream>
#include<vector>
#include<iostream>

write_file::write_file(string nam) {filename=nam;} //cout<<"write obj declared"<<filename; }

void write_file::writep(const vector<double> &vect) {
std::ofstream fout;
fout.open(filename);
for(size_t i=0; i<vect.size(); i++)
{fout<<vect[i]<<std::endl;
}
//cout<<"\n Written successfully";
}

