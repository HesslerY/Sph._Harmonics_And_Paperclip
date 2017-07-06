#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(){
	int numElements = 2664;
	ofstream myfile;
	myfile.open("RSquared.1.txt");
	double theta[numElements];
	double phi[numElements];
	double r[numElements];
	double sphHar[numElements];
	double Pi = M_PI;
	double R2=0.0;
	
// Obtain the r (Gain or Phase) data
	cout << "Enter r: "<< endl;
	for(int i=0; i< numElements; i++){
		cin >> r[i];
	}
	
// Populate the theta array
	for (int j = 0; j < 72; j++){
		for (int k=0; k<37; k++){
			theta[j*37+k]=5*k*3.14159/180;
		}
	}

// Populate the phi array
	for (int l = 0; l < 72; l++){
		for (int m=0; m < 37; m++){
			phi[l*37+m]=5*l*3.14159/180.;
		}
	}
	
// Populate harmonic data points
	for (int p = 0; p < numElements; p++){
		sphHar[p]= 3.14007/(2*sqrt(Pi)) + 
		           (2.026211*pow(10,-6))*cos(theta[p]/ 2)*sqrt(3/Pi)- 
		           0.525095*((3*pow(cos(theta[p]), 2)- 1)/4)*sqrt(5/Pi) + 
				   (7.7912*pow(10, -7))*(1/4)*sqrt(7/Pi)*(5*pow(cos(theta[p]),3)- 3*cos(theta[p]))-
				   0.155673*(3/16)*sqrt(1/Pi)*(35*pow(cos(theta[p]),4) - 30*pow(cos(theta[p]),2)+3)+
				   (7.38295*pow(10,-7))*(sqrt(11/Pi)/16)*(15*cos(theta[p]) - 70*pow(cos(theta[p]),3) + 
  					63*pow(cos(theta[p]),5)) 
		          -0.0954773*(sqrt(13/Pi)/32) *(-5 + 105*pow(cos(theta[p]),2)-315*pow(cos(theta[p]),4) + 231*pow(cos(theta[p]),6))+
				  (-6.50591*pow(10, -7))*(sqrt(15/Pi)/32)*(-35*cos(theta[p])+ pow(315*cos(theta[p]),3) -pow( 693*cos(theta[p]),5) + 429*pow(cos(theta[p]),7))-
				  0.11424*(sqrt(17/Pi)/256)*(35 - 1260*pow(cos(theta[p]),2) + 6930*pow(cos(theta[p]),4) - 12012*pow(cos(theta[p]),6) + 6435*pow((cos(theta[p])),8));
	}

// Calculate R^2
	for(int n=0; n < numElements; n++){
		R2 += pow((r[n]-sphHar[n]),2);
	}
	R2 = R2/numElements;
	
// Close file
	myfile << "Delta R Squared / N: "<< R2 << endl;
	myfile.close();
	return 0;
}
