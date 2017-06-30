#include <iostream>
using namespace std;

int main(){
	int numElements = 0;
	cout << "Number of elements?"<< endl;
	cin >> numElements;
	
	double theta[numElements];
	double phi[numElements];
	double r[numElements];
	
	cout << "Enter theta: "<< endl;

	for(int i=0; i< numElements; i++){
		cin >> theta[i];
	}
	cout << "Enter phi: "<< endl;

	for(int i=0; i< numElements; i++){
		cin >> phi[i];
	}cout << "Enter r: "<< endl;

	for(int i=0; i< numElements; i++){
		cin >> r[i];
	}
	for(int i=0; i< numElements; i++){
		cout <<"{{"<<theta[i]<< ", "<<phi[i]<<"},"<<r[i]<<"},";
	}
	
}
