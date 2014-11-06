#include <iostream>
using namespace std;

int add(int s,int y) {
	return s+y;
}

void Derpy() {
	int p=7;
	cout << p;

}

int Boolcheck(bool bval) {
	if (bval)
		cout << "bval was true" << endl;
	else
		cout << "bval was false" << endl;

}

int forloops(int istart,int iend) {
	for (int iii=istart; iii<iend;iii++)
	return iend+1;

}


int main()
{
	// using namespace std
        cout << "Enter a number: ";
	int x;
	cin >> x;
	cout << "You Entered " << x << endl;
	Derpy(); 
	cout << endl;
	cout << "This is over" << endl;
	cout << add(5,7) << endl;
	Boolcheck(true);
	cout << endl;
	cout << "Okay, now its really over" << endl;
	cout<< forloops(4,9) << endl;
	cout<< endl;
	for (int j=1; j<20;j++)
		cout << j << " ";
	cout<< endl;
	int var=20;
	int *ip;
	ip = &var;
	cout << "Address: " << endl;
	cout << ip << endl;
	cout << "Value: " << endl;
	cout << *ip << endl;
	int rays[4]={1,2,3,4};
	int *r[4];
	int k;
	for (k=0;k<4;k++)
		r[k] = &rays[k];
		cout<< r[k];
	cout << endl;
	int i;
	for (i=0;i<4;i++)
		cout<< rays[i] << " ";
	cout << endl;
	int mat[2][2];
	int m,n;
	for (n=0;n<2;n++)
		for (m=0;m<2;m++)
			mat[m][n] = m*n;
			cout << mat[m][n] << "p ";
	cout << endl;	
	cout<< mat[1][1] << endl;
}


