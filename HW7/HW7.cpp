#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <istream>
#include <list>
#include <string>
#include <math.h>
#include <iterator>


class Matrix {
	public:
	int col, row;
	int datapts;
	std::vector< std::vector<int> > Data;
	std::vector<int> rowstart;
	std::vector<int> rowend;
	void inputdata(char*);
	int values() {return row*col;}
	std::vector<double> multiply(std::vector<double>);
		
};


int printnodes(std::vector<double> x, int nodes) {
	std::vector<double> xtemp(x.size());
	if(nodes<x.size()) {
	for(int i=0;i<x.size();i++) {
		xtemp[i]=x[i];
	}
	for(int j=1;j<=nodes;j++) {
		std::vector<double>::iterator p = std::max_element(xtemp.begin(),xtemp.end());
		if(p!=xtemp.end()) {
			std::vector<double>::iterator its;
			int val = distance(xtemp.begin(),p);
			std::cout<<"Rank: "<<j<<", " << *p << ", "<< val+1 <<std::endl;
			xtemp[val]=0;
		}
	}
	}
	else{std::cout << "Node request is too large"<<std::endl;}
	return 0;
}


std::vector<double> Matrix::multiply(std::vector<double> x) {
	std::vector<double> Ax(row);
	Ax.assign(row,0);
	#pragma omp parallel for private(Ax)
	for(int j=0;j<row;j++){
		if(rowstart[j]!=Data.size()){
			for(int i=rowstart[j];i<rowend[j];i++){
				Ax[j]+=Data[i][2]*x[Data[i][1]];
			}
		}
	}
	return Ax;

}

std::vector<double> normit(std::vector<double> x) {
	int h=x.size();
	double norm = 0;
	for(int i=0;i<h;i++) {
		norm+=x[i]*x[i];
	}
	norm=sqrt(norm);
	for(int i=0;i<h;i++) {
		x[i]=x[i]/norm;
	}
	return x;

}

std::vector<double> randvector(int n,int flag) {
	//flag of 1 will give a vector of all ones;
	//anything else will give a "random" vector;
	std::vector<double> randvec(n);
	if(flag==1){
	for(int i=0;i<n;i++) {
		randvec[i]=1;
	}

	}
	else{
	for(int i=0;i<n;i++) {
		randvec[i]=(i+i/13+.45)/n;
	}
	std::srand(std::time(0));
	std::random_shuffle(randvec.begin(),randvec.end());

	}
return randvec;
}


void Matrix::inputdata(char* fileloc) {

	std::ifstream file;
	file.open(fileloc);
	std::string line;
	col=0;
	row=0;
	int iters=0;
	std::vector<std::string> vals;
	if(file.is_open()) {
		while(std::getline(file,line)) {
			Data.push_back(std::vector<int>());
			std::stringstream stream(line);
			int val1, val2, val3;
			stream >> val1 >> val2 >> val3;
			(Data.back()).push_back(val1); 
			(Data.back()).push_back(val2-1); 
			(Data.back()).push_back(val3);	
			if(val1>row){
			
				row=val1;
			}
			if(val2>col) {
			
				col=val2;
			}
			iters++;
		}
	file.close();
	}
	else{ std::cout<< "File name error - didn't open" << std::endl;}
	datapts=iters;
	std::cout<< datapts<<std::endl;
	std::vector< std::vector<int> > tempdata(3,std::vector<int>(Data.size()));

	for(int i=0;i<datapts;i++) {
		for(int j=0;j<3;j++) {
			tempdata[j][i]=Data[i][j];
		}
	}
	std::cout << "Done Loading Data, start row sort" << std::endl;
	        std::vector<int>::iterator it; 
		rowstart.assign(row,0);
		rowend.assign(row,0);
//
			#pragma omp parallel for private(it)
			for(int i=row-1;i>=0;i--) {
				it = std::lower_bound(tempdata[0].begin(),tempdata[0].end(),i+1);
				rowstart[i] = it - tempdata[0].begin();    
				rowend[i] = rowstart[i];
			}    
					     
		rowend.erase(rowend.begin());
		rowend.push_back(Data.size());
		std::cout << "Finished building row sort" << std::endl;

}


int main(int argc,char* argv[]) {
	
	int T=atoi(argv[1]);
	char* fileloc = argv[3];
	int n_thread=atoi(argv[2]);
	int nodes=atoi(argv[4]);

	Matrix Tester;
	Tester.inputdata(fileloc);

	std::cout<< Tester.row << " " << Tester.col << std::endl;
	omp_set_num_threads(n_thread);

	//Problem 2
	//
	//
	double avgtime=0;
	int iter=1;
	std::vector<double> rando=randvector(Tester.row,0);
	for(int i=1;i<=iter;i++){

		double start=omp_get_wtime();
		rando=Tester.multiply(rando);
		double final=omp_get_wtime() - start;
		avgtime+=final;
		std::cout << i << std::endl;
	}
	avgtime=avgtime/iter;
	std::cout<< "Average Multiply time: "<< avgtime<< std::endl;
	//Problem 3

	double starteig = omp_get_wtime();
	std::vector<double> x=randvector(Tester.row,1);
	x=normit(x);
	for(int i=1;i<=T;i++) {
		x=Tester.multiply(x);
		x=normit(x);
	}

	std::vector<double> xhold(x.size());
	for(int j=0;j<x.size();j++) {
		xhold[j]=x[j];	
	}
	x=Tester.multiply(x);
	//x=normit(x);
	double lam=0;
	for(int i=0;i<x.size();i++){
		lam+=xhold[i]*x[i];
	}	
	double endeig=omp_get_wtime() - starteig;
	std::cout<<"Max Eigenvalue: "<< lam << std::endl;
	std::cout << "Eigenvalue solve time: " << endeig <<std::endl;
	int s = printnodes(xhold,nodes);	

	return 0;
}

