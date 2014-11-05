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

double prederror(std::vector<double> w, std::vector<int> y, std::vector<std::vector<double> > testdata, std::vector<std::vector<int> > feattrain) {
	int count=0;
	for(int i=0;i<y.size();i++) {
		int h=testdata[i].size();
		double tempy=0;
		for(int j=0;j<h;j++) {
			tempy+=w[feattrain[i][j]]*testdata[i][j];
				
		}
		if(y[i]*tempy>0){
			count+=1;
		}
	}

	double error=((double)count/(double)y.size())*100;
	return error;

}

double primval(double C, std::vector<int> y, std::vector<double> w, std::vector<std::vector<double>> testdata, std::vector<std::vector<int>> feattrain) {

	double ww=0;
	double tornado=0;
	for(int i=0;i<w.size();i++) {
		ww+=w[i]*w[i];
	}
	for(int j=0;j<y.size();j++) {

		double temp=0;
		double wx=0;
		double h=testdata[j].size();
		for(int k=0;k<h;k++) {

			wx+=w[feattrain[j][k]]*testdata[j][k];
		}

		double temp1=0;
		temp1=std::max((double)0,(1-y[j]*wx));
		tornado+=(temp1*temp1);
	}
	double fval=(double)ww/2 + (double)C*tornado;
return fval;
}

double wnorm(std::vector<int> y, std::vector<double> alpha, std::vector<std::vector<double>> traindata, int maxn, std::vector<double> w, std::vector<std::vector<int>> feattrain) {
	std::vector<double> wtemp(w.size());
	std::fill(wtemp.begin(),wtemp.end(),0);
	for(int i=0;i<y.size();i++) {
		int h=feattrain[i].size();
		for(int j=0;j<h;j++){
			wtemp[feattrain[i][j]]+=alpha[i]*y[i]*traindata[i][j];
		}	
	}
	double norm=0;
	for(int k=0;k<w.size();k++) {
		wtemp[k] -= w[k];
	}
	for(int k=0;k<w.size();k++) {
		norm+= wtemp[k]*wtemp[k];
}
return norm;
}


double dualval(double C, std::vector<double> alpha, std::vector<double> w) {
	double aa=0;
	double asum=0;
	double aqa=0;
	double fval=0;
	for(int i=0;i<alpha.size();i++) {
		aa=aa + alpha[i]*alpha[i];
		asum=asum +alpha[i];
	}
	for(int j=0;j<w.size();j++){
		aqa+=w[j]*w[j];
	}	
	fval =aqa/2 + aa/(4*C) - asum;
return fval;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
	 std::vector<std::string> elems;
	 split(s, delim, elems);
	 return elems;
}

int main(int argc, char* argv[]) {

	if(argc!=5) {

		std::cout << "You have the wrong argument count" << std::endl;
	}
	
	double C=atof(argv[1]);
	int nr_threads=atoi(argv[2]);
	char*  trainloc = argv[3];
	char*  testloc = argv[4];


	omp_set_num_threads(nr_threads);


	std::vector<int> y;
	std::vector< std::vector<double> > traindata;
	std::vector< std::vector<int> > feattrain;
	std::ifstream file;
	file.open(trainloc);
	std::string line;
	std::vector<std::string> vals;
	std::vector<std::string> vals2;
	std::vector<std::string> feat_val;
	std::vector<int> yt;
	std::vector< std::vector<double> > testdata;
	std::vector< std::vector<int> > feattest;
	std::vector<std::string> feat_valt;
	if(file.is_open()) {
	
		while(std::getline(file,line)) {
			feattrain.push_back(std::vector<int>());
			traindata.push_back(std::vector<double>());

			vals= split(line, ' ');
			y.push_back(std::stoi(vals[0]));
			for(int i=1;i<vals.size();i++) {
				feat_val = split(vals[i],':');
				(feattrain.back()).push_back(std::stoi(feat_val[0])-1);
				(traindata.back()).push_back(std::stof(feat_val[1]));

			}
				
		}
	file.close();
	std::ifstream file2;
	file2.open(testloc);
	std::string line2;
		while(std::getline(file2,line2)) {
			feattest.push_back(std::vector<int>());
			testdata.push_back(std::vector<double>());

			vals2= split(line2, ' ');
			yt.push_back(std::stoi(vals2[0]));
			for(int j=1;j<vals2.size();j++) {
				feat_valt = split(vals2[j],':');
				(feattest.back()).push_back(std::stoi(feat_valt[0])-1);
				(testdata.back()).push_back(std::stof(feat_valt[1]));

			}
				
		}
	file2.close();
	}
	else {
		std::cout << "File didn't open" << std::endl;
	}


	std::cout << "Data has been loaded" << std::endl;

	int sizet=yt.size();
	int size=y.size();
	std::vector<int> permutes(size);
	std::vector<int> used(size);

#pragma omp parallel for
	for(int i=0;i<size;i++) {
		permutes[i]=i;
	}


// permuting the vector and finding the max feature number	
	int maxn=0;
#pragma omp parallel for
	for(int t=0;t<size;t++) {
		int h=traindata[t].size();
#pragma omp parallel for
		for(int ins=0;ins<h;ins++) {
			if(feattrain[t][ins]>maxn){
				maxn=feattrain[t][ins];
			}	
		}
	}
	std::cout << "Max feature count: " << maxn << std::endl;

	std::vector<double> alpha(size), w(maxn+1);
	std::fill(alpha.begin(),alpha.end(),0);
	std::fill(w.begin(),w.end(),0);
	double delt=0;
	for(int outer=1;outer<=20;outer++){
		std::srand(std::time(0));
		std::random_shuffle(permutes.begin(),permutes.end());
		double itertime1 = omp_get_wtime();
		#pragma omp parallel for private(delt) shared(alpha,w,traindata,feattrain) schedule(dynamic) 	
		for(int iter=0;iter<size;iter++) {


			int i=permutes[iter];
			int tempsize=traindata[i].size();
			// x w inner product
			double inprodxw=0;
			double inprodxx=0;
			for(int xal=0;xal<tempsize;xal++){
				inprodxw+=traindata[i][xal]*w[feattrain[i][xal]];
				inprodxx+=traindata[i][xal]*traindata[i][xal];	
			}
			delt=(1-y[i]*inprodxw-((double)alpha[i]/(2*C)))/(y[i]*y[i]*inprodxx + ((double)1/(2*C)));
			double gz = std::max(delt,-alpha[i]);
			alpha[i]+=gz;

			for(int xal2=0;xal2<tempsize;xal2++){
				#pragma omp atomic
				w[feattrain[i][xal2]] = w[feattrain[i][xal2]] + y[i]*gz*traindata[i][xal2];
			}
}
		
		double itertime2 = omp_get_wtime() - itertime1;
		double norms = wnorm(y,alpha,traindata,maxn,w,feattrain);
		double primvals = primval(C,y,w,traindata,feattrain);
		double dualvals = dualval(C,alpha,w);


		double midaccuracy=prederror(w,yt,testdata,feattest);
		std::cout<< outer << "- "<<  "Accuracy: "<<midaccuracy <<", Wtime: " << itertime2<< ", Primal: "<<primvals<<", W norm: "<<norms<<", Dual: "<< dualvals<<std::endl;
	}
	std::cout<< "Finished" << std::endl;

	double accuracy=prederror(w,yt,testdata,feattest);
	std::cout<< accuracy << std::endl;
return 0;
}

