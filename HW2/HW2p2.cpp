#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <istream>
#include"Eigen3.2.2/Eigen/Dense"
#include"Eigen3.2.2/Eigen/Sparse"
#include"Eigen3.2.2/Eigen/Core"


double RMSEcalc(std::vector< std::vector<int> > traindata, Eigen::MatrixXf U, Eigen::MatrixXf M) {
	Eigen::VectorXf sqvals(traindata[1].size());
#pragma omp parallel for
	for(int p=0;p<traindata[1].size();p++) {
		sqvals(p)=(U.col(traindata[0][p]-1).transpose()*M.col(traindata[1][p]-1) - traindata[2][p])*(U.col(traindata[0][p]-1).transpose()*M.col(traindata[1][p]-1) - traindata[2][p]);	
	}	
	double total;
       	total = sqvals.sum();
	return total;


}



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



int main(int argc,char* argv[]) {
	if(argc!=5) {
		std::cout << "You have the wrong argument count" << std::endl;

	}

	//pulling in the variables
	int rank=atoi(argv[1]);
	double lambda=atof(argv[2]);
	int nr_threads=atoi(argv[3]);
	char* data_dr =argv[4];
	char* met = argv[4];
	char trainloc[70];
	char testloc[70];
	strcpy(testloc,argv[4]);
	strcpy(trainloc,argv[4]);

	omp_set_num_threads(nr_threads);
	
	strcat(met,"meta");
	int m,n,trains,tests;
	std::string testname,trainname;

	//input meta file
	std::ifstream metafile;
	metafile.open(met);
	std::string s_line;
	if (metafile.is_open()) {
		std::cout << "Loading Data from MetaFile" << std::endl;
		std::string temp;
		int tempi=0;
		std::getline(metafile,s_line);//getting first line info
		std::stringstream stream(s_line);		
		stream >> tempi;
		m=tempi;
		stream >> tempi;
		n=tempi;	

		std::getline(metafile,s_line);//getting second line
		std::stringstream stream2(s_line);
		tempi=0;
		stream2 >> tempi;
		trains = tempi;
		stream2 >> temp;
		trainname = temp;
	
		std::getline(metafile,s_line);//getting third line
		std::stringstream stream3(s_line);
		tempi=0;
		stream3 >> tempi;
		tests = tempi;
		stream3 >> temp;
		testname = temp;
	
		std:: cout << "Finished Loading Data" << std::endl;	

		metafile.close();
	}
	else {
		std::cout << "file extension incorrect" << std::endl;
	}
	
	strcat(trainloc,trainname.c_str()); //setting location of train file
	strcat(testloc,testname.c_str()); //setting location of test file
	
	std::vector< std::vector<int> > traindata;
	traindata = getdata(trainloc,trains);	
	
	std::vector< std::vector<int> > testdata;
	testdata = getdata(testloc,tests);
	

	std::vector< std::vector<int> > traintemp(trains,std::vector<int>(3));
	std::vector< std::vector<int> > traindatah(3,std::vector<int>(trains));
	
	#pragma omp parallel for
	for(int i=0;i<trains;i++) {
		for(int j=0; j<3;j++) {
			traintemp[i][j]=traindata[j][i];
		}

	}
	//stable_sort
	std::stable_sort(traintemp.begin(),traintemp.end());
	#pragma omp parallel for	
	for(int i=0;i<trains;i++) {
		for(int j=0; j<3;j++) {
			traindatah[j][i]=traintemp[i][j];
		}

	}	
	

	//Find start and end indicies for nonzero entries	
	
	std::vector<int>::iterator it;
	std::vector<int> colstart(n);
	std::vector<int> colend(n);

	#pragma omp parallel for private(it)
	for(int i=n-1;i>=0;i--) {
		it = std::lower_bound(traindata[1].begin(),traindata[1].end(),i+1);
		colstart[i] = it - traindata[1].begin();	
		colend[i] = colstart[i];
	}	
	
	colend.erase(colend.begin());
	colend.push_back(trains);

	std::vector<int>::iterator itr;
	std::vector<int> rowstart(m);
	std::vector<int> rowend(m);

	#pragma omp parallel for private (itr)
	for(int j=m-1;j>=0;j--) {
		itr=std::lower_bound(traindatah[0].begin(),traindatah[0].end(),j+1);
		rowstart[j] = itr - traindatah[0].begin();
		rowend[j] = rowstart[j];
	}
	
	rowend.erase(rowend.begin());
	rowend.push_back(trains);

	std::cout << "Finished evaluating input train data" << std::endl;

	//Initializing U M
	
	Eigen::MatrixXf U = Eigen::MatrixXf::Random(rank,m);
	U=U.cwiseAbs();
	Eigen::MatrixXf M = Eigen::MatrixXf::Zero(rank,n);
	Eigen::MatrixXf Iden = Eigen::MatrixXf::Identity(rank,rank);
	Eigen::MatrixXf Input1 = Eigen::MatrixXf::Zero(rank,rank);
	Eigen::VectorXf Input2(rank);
	Iden = lambda*Iden;
	
	#pragma omp parallel for
	for(int ii=0;ii<m;ii++) {
		if (rowstart[ii]==trains) {
			U.col(ii).setZero();
		}
	}	

	//Beginning of iterator loop

	double RMSETestI = sqrt(RMSEcalc(testdata,U,M)/tests);
	std::cout << "Test RMSE Initial: " << RMSETestI << std::endl;
	
	double starttime = omp_get_wtime();

	for(int UMiter=1; UMiter<=10;UMiter++){

		double startiter = omp_get_wtime();
		std::vector<int> recolu;
		std::vector<int> rvals;
		std::vector<int> rerowm;
		Eigen::MatrixXf Mcol;
		Eigen::MatrixXf Ucol;
		Eigen::VectorXf Rates;

		


		//M iteration
		#pragma omp parallel for private(recolu,rvals,Rates,Ucol,Input1,Input2)	
		for(int jj=0; jj<n; jj++) {
			if(colstart[jj]!=trains) {
				int sizechange=colend[jj]-colstart[jj];
				recolu.resize(sizechange);		
				rvals.resize(sizechange);
				Rates.resize(sizechange);
				Ucol.resize(rank,sizechange);
				std::copy(traindata[0].begin()+colstart[jj],traindata[0].begin()+colend[jj],recolu.begin());
				std::copy(traindata[2].begin()+colstart[jj],traindata[2].begin()+colend[jj],rvals.begin());
				for(int z=0;z<sizechange;z++) {
					Rates(z) = rvals[z];
					Ucol.col(z) = U.col(recolu[z]-1); 
					
				}
				Input1 = (Ucol*(Ucol.transpose()) + Iden);	
				Input2 = Ucol*Rates;
				M.col(jj) = Input1.ldlt().solve(Input2);	
				
			}

		}

		//U loop :)
		#pragma omp parallel for private(rerowm,rvals,Rates,Mcol,Input1,Input2)
		for(int kk=0; kk<m; kk++) {
			if(rowstart[kk]!=trains) {
				int sizechangem=rowend[kk]-rowstart[kk];
				rerowm.resize(sizechangem);		
				rvals.resize(sizechangem);
				Rates.resize(sizechangem);
				Mcol.resize(rank,sizechangem);
				std::copy(traindatah[1].begin()+rowstart[kk],traindatah[1].begin()+rowend[kk],rerowm.begin());
				std::copy(traindatah[2].begin()+rowstart[kk],traindatah[2].begin()+rowend[kk],rvals.begin());
				for(int q=0;q<sizechangem;q++) {
					Rates(q) = rvals[q];
					Mcol.col(q) = M.col(rerowm[q]-1); 
				}	
				

				Input1 = (Mcol*(Mcol.transpose()) + Iden);	
				Input2 = Mcol*Rates;
				U.col(kk) = Input1.ldlt().solve(Input2);	
		
			}

		}
	
		double timeriter = omp_get_wtime() - startiter;
		std::cout << "Iteration: " << UMiter << ", Wall Time: " << timeriter << ", ";

		//double RMSEinter = sqrt(RMSEcalc(traindata,U,M)/trains);
		//std::cout <<"RMSE: "<< RMSEinter << std::endl;
		
		double RMSETest = sqrt(RMSEcalc(testdata,U,M)/tests);
		std::cout << "Test RMSE: " << RMSETest << std::endl;
	}
	
	//output the values needed for testing

	std::cout << "Number of threads: " << nr_threads << std::endl;

	double totaltime = omp_get_wtime() - starttime;
	std::cout << "Total Run time: " << totaltime << std::endl;

	//double RMSETrain;
       	//RMSETrain= sqrt(RMSEcalc(traindata,U,M)/trains);
	//std::cout << "Train RMSE: " << RMSETrain << std::endl;

	double RMSETest = sqrt(RMSEcalc(testdata,U,M)/tests);
	std::cout << "Test RMSE: " << RMSETest << std::endl;


	return 0;
}
