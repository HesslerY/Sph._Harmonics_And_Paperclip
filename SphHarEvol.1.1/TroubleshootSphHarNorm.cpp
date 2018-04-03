
// Written by: Suren Gourapura
// Date: 10/28/17


#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <vector>
using namespace std;

double Pi = 3.141592654;
int SphHarMAX = 13;
int PopMAX = 100;


// Code modified from https://codereview.stackexchange.com/questions/110793/insertion-sort-in-c
void insertionSort(double array[], int length){
    int i,j;
    for (i = 1; i < length; i++) {
        double temp = array[i];
        for (j = i; j > 0 && array[j - 1] < temp; j--) {
            array[j] = array[j - 1];
        }
        array[j] = temp;
    }
}


double integrateSpecies(double species[], int SphHarMAX, double Theta, double Spread){
	
	double Min, Max;
// Make sure that the maximum is less than or equal to 180 degrees and the minimum is greater than or equal to 0
	if (Theta + Spread/2.0 > 180){Max = 180;}
	else{Max = Theta + Spread/2.0;}
			
	if (Theta - Spread/2.0 < 0){Min = 0;}
	else{Min = Theta - Spread/2.0;}

// We need to only integrate over the specified number of spherical harmonics, not all of them. So, we store which spherical harmonics to turn on as an additional binary coefficient.
	int Sph[13];
	for (int i = 0; i < 13; i++){
		if (i < SphHarMAX){
			Sph[i]=1;
		}
		else{
			Sph[i]=0;
		}
	}

	double Sum = 0;			
// Sum the spherical harmonic values by eighth-degree increments from the Min to Max theta values
	for (double degree = Min * Pi / 180; degree <= Max *Pi / 180; degree = degree + Pi/1440){
		Sum += pow(
//Y00		
			Sph[0]*species[0]*(1/2.0)*(1/sqrt(Pi)) 
//Y01
			+Sph[1]*species[1]*(1/2.0)*sqrt(3/Pi)*cos(degree)
//Y02
			+Sph[2]*species[2]*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(degree), 2) - 1)
//Y03
			+Sph[3]*species[3]*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(degree),3) - 3*cos(degree))
//Y04
			+Sph[4]*species[4]*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(degree),4) - 30*pow(cos(degree),2) + 3)
//Y05
			+Sph[5]*species[5]*(1/16.0)*sqrt(11/Pi)*(63*pow(cos(degree),5) - 70*pow(cos(degree),3) + 15*cos(degree))
//Y06
			+Sph[6]*species[6]*(1/32.0)*sqrt(13/Pi)*(231*pow(cos(degree),6) - 315*pow(cos(degree),4) + 105*pow(cos(degree),2) - 5 )
//Y07					
			+Sph[7]*species[7]*(1/32.0)*sqrt(15/Pi)*(429*pow(cos(degree),7) - 693*pow(cos(degree),5) + 315*pow(cos(degree),3) - 35*cos(degree))
//Y08					
			+Sph[8]*species[8]*(1/256.0)*sqrt(17/Pi)*(6435*pow((cos(degree)),8) - 12012*pow(cos(degree),6) + 6930*pow(cos(degree),4) - 1260*pow(cos(degree),2) + 35)
//Y09
			+Sph[9]*species[9]*(1/256.0)*sqrt(19/Pi)*(12155*pow((cos(degree)),9) - 25740*pow(cos(degree),7) + 18018*pow(cos(degree),5) - 4620*pow(cos(degree),3) + 315*cos(degree))
//Y10					
			+Sph[10]*species[10]*(1/512.0)*sqrt(21/Pi)*(46189*pow(cos(degree),10) - 109395*pow((cos(degree)),8) + 90090*pow(cos(degree),6) - 30030*pow(cos(degree),4) + 3465*pow(cos(degree),2) - 63)		
//Y011
			+Sph[11]*species[11]*(1/512.0)*sqrt(23/Pi)*(88179*pow(cos(degree),11) - 230945*pow(cos(degree),9) + 218790*pow((cos(degree)),7) - 90090*pow(cos(degree),5) + 15015*pow(cos(degree),3) - 693*cos(degree))
//Y012
			+Sph[12]*species[12]*(5/2048.0)*sqrt(1/Pi)*(676039*pow(cos(degree),12) - 1939938*pow(cos(degree),10) + 2078505*pow((cos(degree)),8) - 1021020*pow(cos(degree),6) + 225225*pow(cos(degree),4) - 18018*pow(cos(degree),2) + 231)
			, 2)*sin(degree);
	}
// Finally, we give the summed values of the spherical harmonic from theta min to theta max as the test score (a rough kind of integration by .5 degree increments)			
	return 2*Pi*(0.125*Pi/180.0)*Sum;	
}


int main()
{
	double pop[PopMAX][SphHarMAX];
	double Theta = 90; 
	double Spread = 20;
	
// Generate a test population
	for(int i =0; i < PopMAX; i++){
		for(int j =0; j < SphHarMAX; j++){
			pop[i][j] =2.4 ;
		}
	}
	cout << "Initial pop" << endl;
	for(int i =0; i < PopMAX; i++){
		for(int j =0; j < SphHarMAX; j++){
			cout << pop[i][j] << " ";
		}
		cout << endl;
	}	
	
// Normalize the traditional method	
	double integral = 0;	
	for(int i =0; i < PopMAX; i++){
		integral = integrateSpecies(pop[i], SphHarMAX, 90, 180);
		for(int j =0; j < SphHarMAX; j++){
			pop[i][j] = pop[i][j] / sqrt(integral);
		}	
	}	
	
	cout << "Traditional Norm pop" << endl;
	for(int i =0; i < PopMAX; i++){
		for(int j =0; j < SphHarMAX; j++){
			cout << pop[i][j] << " ";
		}
		cout << endl;
	}	

// Generate the test population again
	for(int i =0; i < PopMAX; i++){
		for(int j =0; j < SphHarMAX; j++){
			pop[i][j] = 2.4;
		}
	}

// Normalizing the mathy method		
	integral = 0;
	for(int i =0; i < PopMAX; i++){
		for(int j =0; j < SphHarMAX; j++){
			integral += pow(pop[i][j], 2);
		}
		for(int j =0; j < SphHarMAX; j++){
			pop[i][j] = pop[i][j] / sqrt(integral);
		}	
		integral = 0;
	}	

	cout << "Mathy Norm pop" << endl;
	for(int i =0; i < PopMAX; i++){
		for(int j =0; j < SphHarMAX; j++){
			cout << pop[i][j] << " ";
		}
		cout << endl;
	}	

return 0;
}

