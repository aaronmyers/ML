#include <iostream>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "eigen3.2.2/Eigen/Dense"
#include "eigen3.2.2/Eigen/Core"

double calc_RMSE(Eigen::MatrixXf U, Eigen::MatrixXf M, std::vector< std::vector<int> > R){
	int R_size = R[0].size();
	int ii,jj;
	double estval,trueval;
	double RMSE = 0;
	Eigen::VectorXf ips(R_size);

	for(int k=0 ; k<R_size;k++){
		ii = R[0][k]-1;
		jj = R[1][k]-1;
		estval = ((U.col(ii)).transpose())*M.col(jj);
		trueval = R[2][k];
	       	ips(k) = (estval - trueval)*(estval - trueval);
	}
	RMSE = ips.sum();
	RMSE /= R_size;
	return sqrt(RMSE);
}


std::vector< std::vector<int> > read_data(char * file_loc, int file_size){
	std::ifstream file;
	file.open(file_loc);
	std::vector< std::vector<int> > dataset(3,std::vector<int>(file_size));
	std::string line;
	int count = 0;
	if(file.is_open()){
		while(std::getline(file,line)){
			std::stringstream stream(line);
			stream >> dataset[0][count] >> dataset[1][count] >> dataset[2][count];
			count++;
		}
		file.close();
	}
	else{
		std::cout << "File location incorrect" << std::endl;
	}
	return dataset;
}


int main (int argc, char* argv[]) {
	if (argc!=5){
		std::cout << "You done goofed! Input structure is rank, lambda, nr_threads, data_dir" << std::endl;
	}
	int rank = atoi(argv[1]);
	double lambda = atof(argv[2]);
	int nr_threads = atoi(argv[3]);
	char * metaloc = argv[4];
	std::string str_line, tr_name, te_name;
	char testloc[40];
	char trainloc[40];	
	strcpy(testloc,argv[4]);
	strcpy(trainloc,argv[4]);	
	strcat(metaloc,"meta");
	int m,n,train_size,test_size;
	omp_set_num_threads(nr_threads);

	//read in meta file
	std::ifstream currfile;
	currfile.open(metaloc);
	if(currfile.is_open()){
		int temp = 0;
		std::string str_temp;

		std::getline(currfile,str_line); //first line
		std::stringstream stream(str_line);
		stream >> m >> n;

		std::getline(currfile,str_line); //second lin
		std::stringstream stream2(str_line);
		stream2 >> train_size >> tr_name;

		std::getline(currfile,str_line); //third line
		std::stringstream stream3(str_line);
		stream3 >> test_size >> te_name;

		currfile.close();
	}
	else{
		std::cout << "Meta file not present!" << std::endl;
	}
	currfile.clear();
	strcat(trainloc,tr_name.c_str());
	strcat(testloc,te_name.c_str());
	tr_name.clear();te_name.clear();

	//Read in data	
	std::cout << "Reading training data into file..." << std::endl;
	std::vector< std::vector<int> > colsort;
	colsort = read_data(trainloc,train_size);
	std::vector< std::vector<int> > test_data;
	std::cout << "Reading test data into file..." << std::endl;
	test_data = read_data(testloc,test_size);
	std::cout << "Done reading data" << std::endl;

	//Sort data
	std::cout << "Sorting data..." << std::endl;
	std::vector< std::vector<int> > tmp_rowsort(train_size,std::vector<int>(3));
	std::vector< std::vector<int> > rowsort(3,std::vector<int>(train_size));
	for (int j=0;j<train_size;j++){
		for (int i=0;i<3;i++){
			tmp_rowsort[j][i] = colsort[i][j];
		}
	}

	std::stable_sort(tmp_rowsort.begin(),tmp_rowsort.end());
	
	for(int j=0;j<train_size;j++){
		for(int i=0;i<3;i++){
			rowsort[i][j] = tmp_rowsort[j][i];
		}
	}
	tmp_rowsort.clear();
	std::cout << "Done sorting data" << std::endl;

	//Pull out nonzeros
	std::cout << "Indexing nonzero entires of each column..." << std::endl;
	std::vector<int>::iterator iter;
	std::vector<int> idx_colstart(n);
	std::vector<int> idx_colend(n);

	for(int col=n-1; col>=0; col--){
		iter = std::lower_bound(colsort[1].begin(),colsort[1].end(),col);
		idx_colstart[col] = (iter - colsort[1].begin());
		idx_colend[col] = idx_colstart[col];
		}
	idx_colend.erase(idx_colend.begin());
	idx_colend.push_back(train_size);
	
	std::cout << "Indexing nonzero entires of each row..." << std::endl;
	std::vector<int> idx_rowstart(m);
	std::vector<int> idx_rowend(m);

	for(int row=m-1; row>=0; row--){
		iter = std::lower_bound(rowsort[0].begin(),rowsort[0].end(),row);
		idx_rowstart[row] = (iter - rowsort[0].begin());
		idx_rowend[row] = idx_rowstart[row];
		}
	idx_rowend.erase(idx_rowend.begin());
	idx_rowend.push_back(train_size);
	std::cout << "Done indexing" << std::endl;

	//Run actual iterations
	//Initialize values
	Eigen::MatrixXf identity = Eigen::MatrixXf::Identity(rank,rank);
	identity = identity*lambda; 
	Eigen::MatrixXf U = Eigen::MatrixXf::Constant(rank,m,2);
	//U = U.cwiseAbs();
	Eigen::MatrixXf M = Eigen::MatrixXf::Zero(rank,n);
	Eigen::VectorXf rhs(rank);
	Eigen::MatrixXf LHS = Eigen::MatrixXf::Zero(rank,rank);
	std::vector<int>::iterator mstart = colsort[0].begin();
	std::vector<int>::iterator ustart = rowsort[1].begin();
	std::vector<int>::iterator rmstart = colsort[2].begin();
	std::vector<int>::iterator rustart = rowsort[2].begin();
	
	//Zero appropriate columns of U
	for (int c = 0; c<m;c++){
		if(idx_rowstart[c] == idx_rowend[c]){
			U.col(c).setZero();
		}
	}
	double time1 = omp_get_wtime();
	std::cout << "Beginning U M iterations..." << std::endl;
	for (int it=0; it<10; it++){	
		std::cout << it << std::endl;	
		std::vector<int> nzidx,rc;
		Eigen::VectorXf r;		
		Eigen::MatrixXf U_chunk;
		Eigen::MatrixXf M_chunk;
	
		//M-step
		for (int mcols=0;mcols<n;mcols++){
			if(idx_colstart[mcols] != idx_colend[mcols]){
				int currlength = idx_colend[mcols] - idx_colstart[mcols];
				U_chunk.resize(rank,currlength);
				r.resize(currlength);
				rc.resize(currlength);
				nzidx.resize(currlength);
				std::copy(idx_colstart[mcols]+mstart,idx_colend[mcols]+mstart,nzidx.begin());
				std::copy(idx_colstart[mcols]+rmstart,idx_colend[mcols]+rmstart,rc.begin());

				for(int nz=0; nz<currlength;nz++){
					U_chunk.col(nz) = U.col(nzidx[nz]-1);
					r(nz) = rc[nz];
				}
				rhs = U_chunk*r;
				LHS = (U_chunk*(U_chunk.transpose())) + identity;
				M.col(mcols) = LHS.ldlt().solve(rhs);
			}
		}

		//U-step
		for (int ucols=0;ucols<m;ucols++){
			if(idx_rowstart[ucols] != idx_rowend[ucols]){
				int currlength = idx_rowend[ucols] - idx_rowstart[ucols];
				M_chunk.resize(rank,currlength);
				r.resize(currlength);
				rc.resize(currlength);
				nzidx.resize(currlength);
				std::copy(idx_rowstart[ucols]+ustart,idx_rowend[ucols]+ustart,nzidx.begin());
				std::copy(idx_rowstart[ucols]+rustart,idx_rowend[ucols]+rustart,rc.begin());

				for(int nz=0; nz<currlength;nz++){
					M_chunk.col(nz) = M.col(nzidx[nz]-1);
					r(nz) = rc[nz];
				}

				rhs = M_chunk*r;
				LHS = (M_chunk*(M_chunk.transpose())) + identity;
				U.col(ucols) = LHS.ldlt().solve(rhs);
			}
		}
		std::cout << calc_RMSE(U,M,colsort) << std::endl;
	}

	std::cout << "Done with U M iterations" << std::endl;
	double train_RMSE = calc_RMSE(U,M,colsort);
	std::cout << "Training RMSE is  " << train_RMSE << std::endl;
	double test_RMSE = calc_RMSE(U,M,test_data);
	std::cout << "Test RMSE is  " << test_RMSE << std::endl;
	std::cout << "Time to factor is " <<omp_get_wtime() -time1 << std::endl;
	return 0;
}

