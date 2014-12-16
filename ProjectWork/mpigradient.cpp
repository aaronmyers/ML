#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <iostream>     
#include <sstream>      
#include <stdio.h>
#include <string.h>
#include <cstdio>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#include <math.h>
#include"eigen-eigen-1306d75b4a21/Eigen/Dense"
#include"eigen-eigen-1306d75b4a21/Eigen/Sparse"
#include"eigen-eigen-1306d75b4a21/Eigen/Core"
#include"eigen-eigen-1306d75b4a21/Eigen/SparseQR"
#include"eigen-eigen-1306d75b4a21/Eigen/QR"
#include"eigen-eigen-1306d75b4a21/Eigen/Householder"
#include"eigen-eigen-1306d75b4a21/Eigen/OrderingMethods"


Eigen::VectorXd grad(Eigen::SparseMatrix<double> Aloct, int m, Eigen::VectorXd RHS, double rho) {

	Eigen::VectorXd xtemp(m);
	Eigen::VectorXd xe(m);
	Eigen::VectorXd diff(m);
	Eigen::VectorXd bas2(m);
	Eigen::VectorXd bas(Aloct.rows());
	xe.fill(4);
	diff.fill(3);
	double alpha = .01;
	while(diff.norm()>.01){
		xtemp = xe;
		bas=Aloct*xtemp;
		bas2 = Aloct.transpose()*bas;	
		xe = xtemp - alpha*(bas2 + rho*xtemp + RHS);
		diff=xe-xtemp;
		std::cout<< diff.norm()<<std::endl;
	}
	

return xe;
}

std::vector<int> getcounts(std::vector<std::vector<int> > Data) {

	std::vector<int> places;
	int val=1;
	int count=0;
	places.push_back(0);
	for(int i=0;i<Data[0].size();i++){
		if(Data[1][i]!=val){
			if(Data[1][i]==val+1){
				places.push_back(count);
				val=val+1;
				count++;
			}
			else{
				while(Data[1][i]!=val+1){
				places.push_back(places[val-1]);
				val++;
				}
				places.push_back(count);
				count++;
				
			}
		}	
		else{count++;}
	}
	places[val+1]=count;

return places;
}

typedef Eigen::Triplet<double> T;

std::vector<double> normit(std::vector<double> x) {
	double sums=0;
	for(int i=0;i<x.size();i++){
		sums+=x[i]*x[i];
	}

	sums=sqrt(sums);
	for(int i=0;i<x.size();i++) {
		x[i]=x[i]/sums;
	}
return x;
}

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
            xtemp[val]=-11342;
        }   
    }   
    }   
    else{std::cout << "Node request is too large"<<std::endl;}
    return 0;
}



std::vector<double> updateu(std::vector<double> utemp,std::vector<double> x,std::vector<double>z) {
	std::vector<double> unew(utemp.size());
	for(int i=0;i<utemp.size();i++) {
		unew[i]=utemp[i]+x[i]-z[i];
	}

return unew;
}



std::vector<double> updatez(std::vector<double> fullu, std::vector<double> fullx, double rho, int size, double lam, int m) {
	std::vector<double> outz(m);
	double comp=lam/(rho*(double)size);
	for(int i=0;i<m;i++) {
		double dumx=0;
		double dumu=0;
		double dux=0;
		for(int j=0;j<size;j++) {
			dumx+=fullx[m*j+i]/(double)size;	
			dumu+=fullu[m*j+i]/(double)size;
		}
			double minus=-1.0;
			dux=dumu+dumx;
			if(dux>comp){
				outz[i]=dux-comp;
			}
			else if(dux<minus*comp) {
				outz[i]=dux+comp;
			}
			else{outz[i]=(double)0;}
	}

return outz;
}


std::vector<double> randvector(int n,double c,int m) {                                                                                    
    //flag of 1 will give a vector of all ones;                                                                                     
    //anything else will give a "random" vector;                                                                                    
    std::vector<double> randvec(n);                                                                                                 
    for(int i=0;i<n;i++) {                                                                                                          
        randvec[i]=(double)(1-c)/(double)m;                                                                                                               
    }
return randvec;                                                                                                                     
}     



std::vector< std::vector<int> >  getdata(char* f_loc) {  
	std::ifstream file;
	file.open(f_loc);
	std::string line;
	std::vector< std::vector<int> > matrix(2,std::vector<int>());
	if(file.is_open()) {
	int count=0;       
	int v1;
	int v2;
	int colcount=1;
	int maxcol=0;
	int m=1;
	while (std::getline(file,line)) {
	    std::stringstream streamd(line);
	    streamd >> v1 >> v2;  
	    matrix[0].push_back(v1);
	    matrix[1].push_back(v2);
	    	
	    if(m<v2) {m=v2;if(colcount>maxcol){maxcol=colcount;} colcount=1;}
	    else{colcount++;}
	    count++;                                                              
	                                                                          
	        }                                                                 
	matrix[1].push_back(m);
	matrix[1].push_back(count);
	matrix[1].push_back(maxcol);
	}                                                                         
	else {MPI_Abort(MPI_COMM_WORLD,1);}
	file.close();
	return matrix;                                                            
}

std::vector<Eigen::Triplet<double> > filltrip(int nnz,int m,std::vector<std::vector<int> > Data , int rank,int size, std::vector<int> places,double c) {
	std::vector<T> triplist;
	int start = places[floor(m/size)*rank];
	int end;
	if(rank==size-1){end = nnz-1;}
	else{end=places[floor(m/size)*(rank+1)]-1;}
	int stopit=end-start+1;
	double third=(double)(1-c)/stopit;
	for(int i=0;i<stopit ;i++) {
		triplist.push_back(T(Data[0][start+i]-1,Data[1][start+i]-1-rank*floor(m/size),third));		
	}
	std::vector<std::vector<int> >().swap(Data);
return triplist;
}

int main(int argc,char* argv[]) {

	char* floc = argv[1];
	int nnz,maxcol;
	double c = atof(argv[2]);
	int maxiter = atoi(argv[3]);
	double rho=atof(argv[4]);
	double lam=atof(argv[5]);
	int rank, size, m;
	MPI_Init(&argc,&argv);
	MPI_Status status;
	MPI_Request request;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

// --------------Loading data (A'), All processors -------------------------------



	std::vector<std::vector<int> > Data=getdata(floc);
	maxcol=Data[1].back();
	Data[1].pop_back();
	nnz=Data[1].back();
	Data[1].pop_back();
	m=Data[1].back();
	Data[1].pop_back();

	if(rank==0){
		printf("Loading the Data...\n");
		printf("[%d] nnz= %d\n",rank,nnz);
		printf("[%d] m= %d\n",rank,m);
		printf("[%d]:Data loaded\n",rank);
	}
	MPI_Barrier(MPI_COMM_WORLD);
//----------------Done loading data still need to modify to get A matrix----	



//---------------- Creating partial bs for each processor-----------------------

	Eigen::VectorXd blocal(floor(m/size));
	if(rank<size-1) {
		blocal.fill((1-c)/m);
	}
	else{blocal.resize(m-rank*floor(m/size)); blocal.fill((1-c)/m);}

//-------------Determining size to allocate for Eigen Matrix-----------------------

	double t1, t2;
	t1 = MPI_Wtime();
	std::vector<int> places=getcounts(Data);

	Eigen::SparseMatrix<double> Aloct(m,floor(m/size)); 
	std::vector<Eigen::Triplet<double> > tripit=filltrip(nnz,m,Data,rank,size,places,c);

	
	if(rank<size-1){
		Aloct.reserve(tripit.size());
		Aloct.setFromTriplets(tripit.begin(),tripit.end());	
		Aloct=Aloct.transpose();
	}	
	else{
		Aloct.resize(m,m-rank*floor(m/size));
		Aloct.reserve(tripit.size());
		Aloct.setFromTriplets(tripit.begin(),tripit.end());	
		Aloct=Aloct.transpose();
	}

//-----------Separate A into each processor and make it (I - c Pt) ------------


	t2 = MPI_Wtime() - t1;
	printf("[%d]: time to load data: %f\n",rank,t2);
	printf("[%d]: size of my matrix: %d x %d\n",rank,Aloct.rows(),Aloct.cols()); 
	int sp=blocal.size();
	printf("[%d]: size of b is: %i\n",rank,sp);
// -----------REMOVE DATA FROM MEMORY-----------------------------------------

	MPI_Barrier(MPI_COMM_WORLD);
//--------------Finished, each processor having its own data ------------------

	Aloct.makeCompressed();
	//Eigen::SparseMatrix<double> FAtA(m,m);
	//FAtA=Aloct.transpose()*Aloct;
	//FAtA.makeCompressed();
	//Eigen::MatrixXd AtA;
	//AtA=Eigen::MatrixXd(FAtA);

	//for(int i=0;i<AtA.rows();i++){
	//	AtA.coeffRef(i,i)+=rho;
	//}
	Eigen::VectorXd Atb(m);
	Atb = Aloct.transpose()*blocal;

	MPI_Barrier(MPI_COMM_WORLD);

//----------Ready for linear solver --------------------------------------


//Initialize z(global) and u(local)
	std::vector<double> z(m);
	if(rank==0){
		z=randvector(m,.3,m);
	}	
	
	MPI_Bcast(&z[0],m,MPI_DOUBLE,0,MPI_COMM_WORLD);
	std::vector<double> u=randvector(m,.6*rank,m);
	//std::vector<double> x=randvector(m,.6*rank,m);//Remove once finished


//---Iterate a set number of times (solve min ||Ax-b||) 
	double finalres;
	int iter=0;
	std::vector<double> fullu(1);
	std::vector<double> fullx(1);
	if(rank==0){fullu.resize(m*size); fullx.resize(m*size);}
	double totime1=MPI_Wtime();

	while(iter<maxiter) {
		//Solve for x update x= (A'A + rho*I)\(A'b + rho*z - rho*u)--------------------------	
		//LHS is now denoted AtA
		//update the RHS of the equation above, Eigen solve the equation, convert x to std::vector

		Eigen::Map<Eigen::VectorXd> ze(&z[0],m);
		Eigen::Map<Eigen::VectorXd> ue(&u[0],m);
		Eigen::VectorXd RHS = rho*ue - rho*ze - Atb;
		Eigen::VectorXd xe(m);
		//Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int> > solver;
		//AtA.makeCompressed();
		//solver.compute(AtA);
		//xe = AtA.householderQr().solve(RHS);
		xe=grad(Aloct,m,RHS,rho);
		std::cout<< xe(1) <<std::endl;
		std::vector<double> x(xe.data(), xe.data() + xe.rows());	
		if(rank==0){
		for(int i=0;i<xe.size();i++){
			std::cout<<xe(i)<<std::endl;
		}
		}
		x=normit(x);

		// z update on the mater thread and then send back out to all other threads----------
		MPI_Gather(&x[0],m,MPI_DOUBLE,&fullx[0],m,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&u[0],m,MPI_DOUBLE,&fullu[0],m,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			z=updatez(fullu,fullx,rho,size,lam,m);	

		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		
		//Send z back to all other threads--------------------------------------------------------
		MPI_Bcast(&z[0],m,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		//Update u now that z is updated ------------------------------
		std::vector<double> utemp(m);
		std::copy(u.begin(),u.end(),utemp.begin());
		u=updateu(utemp,x,z);

		iter++;
		MPI_Barrier(MPI_COMM_WORLD);
		//if(iter==maxiter-1){
		//	Eigen::VectorXd res=Aloct*ze-blocal;
		//	double normval=res.norm();
		//	MPI_Reduce(&normval,&finalres,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		//}
	}	


	
	double finaltime=MPI_Wtime() - totime1;

	std::cout<<"Here come the nodes"<<std::endl;

	if(rank==0){int nothing=printnodes(z,5); printf("[%d]: Total time for %d iterations: %f, Res:\n",rank,maxiter,finaltime);}	
//-----------End of linear solver---------------------------------------

	MPI_Finalize();
return 0;
}
