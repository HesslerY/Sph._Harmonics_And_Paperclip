#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(){
	int numElements = 2664;
	ofstream myfile;
	myfile.open("RSquared.Y00-Y00.txt");
	double theta[numElements];
	double phi[numElements];
	double r[numElements];
	double sphHar[numElements];
	double Pi = M_PI;
	double R2=0.0;
	

	
// Populate the theta array
	fill(theta, theta+numElements, 0.0);
	for (int j = 0; j < 72; j++){
		for (int k=0; k<37; k++){
			theta[j*37+k]=5*k*3.14159/180;
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
		sphHar[p]=3.14007*(1/2.0)*(1/sqrt(Pi));// +2.026211*pow(10,-6)*(1/2)*sqrt(3/Pi)*cos(theta[p]) -0.525095*(1/4)*sqrt(5/Pi)*(3*pow(cos(theta[p]), 2)- 1);
//Y03
//					+7.7912*pow(10, -7)*(1/4)*sqrt(7/Pi)*(5*pow(cos(theta[p]),3)- 3*cos(theta[p]))
//Y04
//					-0.155673*(3/16)*sqrt(1/Pi)*(35*pow(cos(theta[p]),4) - 30*pow(cos(theta[p]),2)+3)
//Y05
//					+7.38295*pow(10,-7)*(1/16)*sqrt(11/Pi)*(15*cos(theta[p]) - 70*pow(cos(theta[p]),3)+63*pow(cos(theta[p]),5)) 
//Y06
//					-0.0954773*(1/32)*sqrt(13/Pi)*(-5 + 105*pow(cos(theta[p]),2)-315*pow(cos(theta[p]),4) + 231*pow(cos(theta[p]),6))
//Y07					
//					-6.50591*pow(10, -7)*(1/32)*sqrt(15/Pi)*(-35*cos(theta[p])+ 315*pow(cos(theta[p]),3) -693*pow(cos(theta[p]),5) + 429*pow(cos(theta[p]),7))
//Y08					
//					-0.11424*(1/256)*sqrt(17/Pi)*(35 - 1260*pow(cos(theta[p]),2) + 6930*pow(cos(theta[p]),4) - 12012*pow(cos(theta[p]),6) + 6435*pow((cos(theta[p])),8))
		cout << theta[p]<<":"<<sphHar[p]<<endl;
	}


// Close file
	myfile << "Delta R Squared / N: "<< R2 << endl;
	myfile.close();
	return 0;
}
