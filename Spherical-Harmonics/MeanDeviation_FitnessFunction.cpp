#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(){
	int numElements = 2664;
	ofstream myfile;
	myfile.open("MeanDeviation.txt");
	double theta[numElements];
	double phi[numElements];
	double r[numElements];
	double sphHar[numElements];
	double Pi = M_PI;
	double md=0.0;
	
// Obtain the r (Gain or Phase) data
	cout << "Enter r: "<< endl;
	fill(r, r+numElements, 0.0);
	for(int i=0; i< numElements; i++){
		cin >> r[i];
	}
	
// Populate the theta array
	fill(theta, theta+numElements, 0.0);
	for (int j = 0; j < 72; j++){
		for (int k=0; k<37; k++){
			theta[j*37+k]=5*k*Pi/180.0;
		}
	}

// Populate the phi array
	fill(phi, phi+numElements, 0.0);
	for (int l = 0; l < 72; l++){
		for (int m=0; m < 37; m++){
			phi[l*37+m]=5*l*3.14159/180.0;
		}
	}
	
// Populate harmonic data points
	fill(sphHar, sphHar+numElements, 0.0);
	for (int p = 0; p < numElements; p++){
		sphHar[p]= 
//Y00		
		3.14007*(1/2.0)*(1/sqrt(Pi)) 
//Y01
		+2.026211*pow(10,-6)*(1/2.0)*sqrt(3/Pi)*cos(theta[p])
//Y02
		-0.525095*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(theta[p]), 2)- 1)
//Y03
		+7.7912*pow(10, -7)*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(theta[p]),3)- 3*cos(theta[p]))
//Y04
		-0.155673*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(theta[p]),4) - 30*pow(cos(theta[p]),2)+3)
//Y05
		+7.38295*pow(10,-7)*(1/16.0)*sqrt(11/Pi)*(15*cos(theta[p]) - 70*pow(cos(theta[p]),3)+63*pow(cos(theta[p]),5))
//Y06
		-0.0954773*(1/32.0)*sqrt(13/Pi)*(-5 + 105*pow(cos(theta[p]),2)-315*pow(cos(theta[p]),4) + 231*pow(cos(theta[p]),6))
//Y07					
		-6.50591*pow(10, -7)*(1/32.0)*sqrt(15/Pi)*(-35*cos(theta[p])+ 315*pow(cos(theta[p]),3) -693*pow(cos(theta[p]),5) + 429*pow(cos(theta[p]),7))
//Y08					
		-0.11424*(1/256.0)*sqrt(17/Pi)*(35 - 1260*pow(cos(theta[p]),2) + 6930*pow(cos(theta[p]),4) - 12012*pow(cos(theta[p]),6) + 6435*pow((cos(theta[p])),8))
		;
	}

// Calculate mean deviation
	for(int n=0; n < numElements; n++){
		md += abs(r[n]-sphHar[n]);
	}
	md = md/numElements;
	
// Close file
	myfile << "Mean Deviation: "<< md << endl;
	myfile.close();
	return 0;
}
