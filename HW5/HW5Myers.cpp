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
#include "Galois/Galois.h"
#include "Galois/Graph/Graph.h"

typedef Galois::Graph::LC_CSR_Graph<double, double> Graph;
typedef Graph::GraphNode GNode;
typedef std::pair<int, GNode> Task;

struct P{
	Graph &graph;
	std::vector<double> &x;
	std::vector<double> &y;
	P(Graph &graph, std::vector<double> &x, std::vector<double> &y) : graph(graph), x(x), y(y) {}
	
	void operator() (unsigned int noduh) {
		double sumit=0;
		for(auto edge : graph.out_edges(noduh)) {
			sumit+=(graph.getEdgeData(edge))*x[graph.getEdgeDst(edge)];
		}	
	y[noduh]=sumit;
	}
	
};

std::vector<double> multiply(Graph &gr, std::vector<double> x) {

	std::vector<double> y(x.size());
	std::fill(y.begin(),y.end(),0);
	Galois::do_all(gr.begin(),gr.end(),P(gr,x,y));
return y;
};

int size(Graph& gr) {
	unsigned int sizer=0;
	for(auto node : gr) {
		sizer++;
	}
	return sizer;
};


std::vector<double> scalarmult(std::vector<double> x, double a) {
	for(int i=0;i<x.size();i++) {
		x[i]=x[i]*a;
	}
	return x;
};

std::vector<double> diagmult(std::vector<double> d, std::vector<double> r) {
	std::vector<double> temp(d.size());
	for(int i=0;i<d.size();i++){
		temp[i]=((double)1/d[i])*r[i];
	}
	return temp;
};

double sumvec(std::vector<double> x) {
	double sumit=0;
	for(int i=0;i<x.size();i++) {
		sumit+=x[i];
	}
	return sumit;
};

std::vector<double> vectoradd(std::vector<double> x, std::vector<double> y) {
	std::vector<double> xy(x.size());
	for(int i=0;i<x.size();i++){
		xy[i]=x[i] + y[i];
	}
	return xy;
	
};

void nodes(std::vector<double> r) {
	std::vector<double> rtemp(r.size());
	for(int i=0;i<r.size();i++) {
		rtemp[i]=r[i];
	}
	for(int j=1;j<=10;j++) {
		std::vector<double>::iterator p=std::max_element(rtemp.begin(),rtemp.end());
		std::vector<double>::iterator its;
		int val= distance(rtemp.begin(),p);
		std::cout<<"Rank: "<<j<<", "<< *p<<", "<< val << std::endl;
		rtemp[val]=0;
	
	}
};

int main(int argc, char* argv[]) {

	if(argc!=3){
		std::cout<<"Argument count is wrong: exec threads fileloc"<<std::endl;
	}
	std::clock_t start;
	double duration;
	start=std::clock();
	int n_threads=atoi(argv[1]);
	char* filename = argv[2];
	Graph graph;
	Galois::Graph::readGraph(graph, filename);
	int gsize=size(graph);
	std::vector<double> x(gsize);
	std::fill(x.begin(),x.end(),1);
	int nruns=20;
	std::clock_t starti;
	double durationi;
	starti=std::clock();
	for(int i=1;i<=nruns;i++) {
		std::vector<double> y=multiply(graph,x);	
	}
	durationi=(std::clock() - starti)/(double) CLOCKS_PER_SEC;
	std::cout<< "Threads: "<<n_threads<<", Time: "<< durationi/(double)nruns<<std::endl;
	double alpha=0.15;
	int T=50;
	Galois::setActiveThreads(4);
	std::vector<double> d(gsize);
	for(auto n : graph) {
		unsigned int rowsum=0;
		for(auto edge : graph.out_edges(n)) {
			rowsum+=graph.getEdgeData(edge);

		}	
		d[n]=rowsum;

	}

	std::vector<double> r(gsize);
	std::vector<double> v(gsize);
	std::fill(v.begin(),v.end(),(double)1/(double)gsize);
	std::fill(r.begin(),r.end(),(double)1/(double)gsize);
	for(int i=0; i<=T;i++){
		double rs=sumvec(r);
		std::vector<double> av=scalarmult(v,alpha);
		av=scalarmult(av,rs);
		std::vector<double> adr=diagmult(d,r);
		adr=scalarmult(adr,(1-alpha));
		adr=multiply(graph,adr);
		r=vectoradd(adr,av);
	}

	nodes(r);
return 0;
}

