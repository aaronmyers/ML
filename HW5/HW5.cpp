#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <sstream>
#include <algorithm>i
#include <stdio.h>
#include <string.h>
#include <istream>
#include"Eigen3.2.2/Eigen/Dense"
#include"Eigen3.2.2/Eigen/Sparse"
#include"Eigen3.2.2/Eigen/Core"

std::vector< std::vector<int> >  getdata(char* f_loc,int f_size) {
        std::ifstream file;
        file.open(f_loc);
        std::string line;
        std::vector< std::vector<int> > matrix(3,std::vector<int>(f_size));
        if(file.is_open()) {
        int count=0;
        while (std::getline(file,line)) {
                std::stringstream streamd(line);
                streamd >> matrix[0][count] >> matrix[1][count] >> matrix[2][count];
                count++;    
    
                        }   
        }   
        return matrix;
}


int main() {


return 0;
}
