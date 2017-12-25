// SphHarEvol.1.cpp.cpp : Defines the entry point for the console application.
// Written by: Suren Gourapura
// Date: 10/28/17
// Goal: The goal is to take a given theta angle and spread in theta to evolve an antenna using spherical harmonics that best satisfies it.

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
using namespace std;

/*
void popTester (int PopMAX, double Theta, double Spread, double pop[][], double testScores[]){
// Run the tester function for each species in the population
	for (int i = 0; i <= PopMAX - 1; i++){
		testScores[i] = tester(Theta,Spread,pop[][],i);
	}
}


double tester(double Theta, double Spread,  int i, double pop[i][]){
// For this code, the goal will be to maximize the value of the species from [Theta - Spread/2] to [Theta + Spread/2]
	double Sum = 0;
	double Min = Theta - Spread/2;
	double Max;
// Make sure that the maximum is less than or equal to 180 degrees
	if (Theta + Spread/2 > 180){
		Max = 180;
	}
	else{
		Max = Theta + Spread/2;
	}
// Sum the spherical harmonic values by degree increments from the Min to Max theta values	
	for (double degree = Min; degree <= Max; degree++){
		Sum += sphericalHarmonicVal(degree, pop[i][], i);
	}

	return Sum;
}


double sphericalHarmonicVal(double degree, double pop[i][]){
// Calculate the spherical harmonic value of species i at theta = degree by calculating the linear sum of Y0,0 - Y0,12 with the coefficients described by pop[i][]
	double Value = 
//Y00		
		pop[i][0]*(1/2.0)*(1/sqrt(Pi)) 
//Y01
		+pop[i][1]*pow(10,-6)*(1/2.0)*sqrt(3/Pi)*cos(degree)
//Y02
		+pop[i][2]*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(degree), 2)- 1)
//Y03
		+pop[i][3]*pow(10, -7)*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(degree),3)- 3*cos(degree))
//Y04
		+pop[i][4]*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(degree),4) - 30*pow(cos(degree),2)+3)
//Y05
		+pop[i][5]*pow(10,-7)*(1/16.0)*sqrt(11/Pi)*(15*cos(degree) - 70*pow(cos(degree),3)+63*pow(cos(degree),5))
//Y06
		+pop[i][6]*(1/32.0)*sqrt(13/Pi)*(-5 + 105*pow(cos(degree),2)-315*pow(cos(degree),4) + 231*pow(cos(degree),6))
//Y07					
		+pop[i][7]*pow(10, -7)*(1/32.0)*sqrt(15/Pi)*(-35*cos(degree)+ 315*pow(cos(degree),3) -693*pow(cos(degree),5) + 429*pow(cos(degree),7))
//Y08					
		+pop[i][8]*(1/256.0)*sqrt(17/Pi)*(35 - 1260*pow(cos(degree),2) + 6930*pow(cos(degree),4) - 12012*pow(cos(degree),6) + 6435*pow((cos(degree)),8))
//Y09
		+pop[i][9]*(1/256.0)*sqrt(19/Pi)*(315*cos(degree)- 4620*pow(cos(degree),3) + 18018*pow(cos(degree),5) - 25740*pow(cos(degree),7) + 12155*pow((cos(degree)),9))
//Y10					
		+pop[i][10]*(1/512.0)*sqrt(21/Pi)*(-63 +3465*pow(cos(degree),2) - 30030*pow(cos(degree),4) + 90090*pow(cos(degree),6) -109395*pow((cos(degree)),8)+46189*pow(cos(degree),10))		
//Y011
		+pop[i][11]*(1/512.0)*sqrt(23/Pi)*(-693*pow(cos(degree),1) +15015*pow(cos(degree),3) - 90090*pow(cos(degree),5) +218790*pow((cos(degree)),7)-230945*pow(cos(degree),9)+88179*pow(cos(degree),11))
//Y012
		+pop[i][12]*(1/2048.0)*sqrt(25/Pi)*(231 -18018*pow(cos(degree),2) +225225*pow(cos(degree),4) - 1021020*pow(cos(degree),6) +2078505*pow((cos(degree)),8)-1939938*pow(cos(degree),10)+676039*pow(cos(degree),12))
		;
	
	return Value;
}
*/

int main()
{
	double Pi = 3.14159;
	int SphHarMAX = 13;
	int PopMAX = 100;
	double pop[PopMAX][SphHarMAX];
	double Theta, Spread;
	int Gen = 0;
	srand(time(NULL));
	
// Enter the goal theta and spread
	cout << "Enter the goal theta and spread (in degrees), and generations in the following format: theta spread generations" << endl;
	cin >> Theta;
	cin >> Spread;
	cin >> Gen;
	
// Randomly fill the population. First step is to give each value in each species random number between 0 and 1. Then, we normalize the species values to 1
	for (int i = 0; i <= PopMAX - 1; i++){
		for (int j = 0; j <= SphHarMAX - 1; j++){
			pop[i][j] = ((double)rand() / (double)(RAND_MAX)); // gives a random value between 0 and 1
		}
	}
	
// We add up the values for each species and divide them by the sum to normalize to 1
	for (int i = 0; i <= PopMAX - 1; i++){
		double counter = 0;
		for (int j = 0; j <= SphHarMAX - 1; j++){
			counter += pop[i][j];
		}
		for (int j = 0; j <= SphHarMAX - 1; j++){
			pop[i][j] = pop[i][j] / counter;
		}		
	}

// Now, we evolve. We loop the following code for each generation
	for (int g = 1; g <= Gen; g++){
		
		cout << "Generation: " << g << endl;
		
// The first objective is to calculate how well our current population's species are doing
		double testScores[100];
		
		for (int i = 0; i <= PopMAX - 1; i++){
				
// For this code, the goal will be to maximize the value of the species from [Theta - Spread/2] to [Theta + Spread/2]
			double Sum = 0;
			double Min = Theta - Spread/2.0;
			double Max;
			
// Make sure that the maximum is less than or equal to 180 degrees
			if (Theta + Spread/2.0 > 180){Max = 180;}
			else{Max = Theta + Spread/2.0;}
			
// Sum the spherical harmonic values by degree increments from the Min to Max theta values	
			for (double degree = Min * Pi / 180; degree <= Max *Pi / 180; degree = degree + Pi/180){
				Sum += 
//Y00		
					pop[i][0]*(1/2.0)*(1/sqrt(Pi)) 
//Y01
					+pop[i][1]*pow(10,-6)*(1/2.0)*sqrt(3/Pi)*cos(degree)
//Y02
					+pop[i][2]*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(degree), 2)- 1)
//Y03
					+pop[i][3]*pow(10, -7)*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(degree),3)- 3*cos(degree))
//Y04
					+pop[i][4]*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(degree),4) - 30*pow(cos(degree),2)+3)
//Y05
					+pop[i][5]*pow(10,-7)*(1/16.0)*sqrt(11/Pi)*(15*cos(degree) - 70*pow(cos(degree),3)+63*pow(cos(degree),5))
//Y06
					+pop[i][6]*(1/32.0)*sqrt(13/Pi)*(-5 + 105*pow(cos(degree),2)-315*pow(cos(degree),4) + 231*pow(cos(degree),6))
//Y07					
					+pop[i][7]*pow(10, -7)*(1/32.0)*sqrt(15/Pi)*(-35*cos(degree)+ 315*pow(cos(degree),3) -693*pow(cos(degree),5) + 429*pow(cos(degree),7))
//Y08					
					+pop[i][8]*(1/256.0)*sqrt(17/Pi)*(35 - 1260*pow(cos(degree),2) + 6930*pow(cos(degree),4) - 12012*pow(cos(degree),6) + 6435*pow((cos(degree)),8))
//Y09
					+pop[i][9]*(1/256.0)*sqrt(19/Pi)*(315*cos(degree)- 4620*pow(cos(degree),3) + 18018*pow(cos(degree),5) - 25740*pow(cos(degree),7) + 12155*pow((cos(degree)),9))
//Y10					
					+pop[i][10]*(1/512.0)*sqrt(21/Pi)*(-63 +3465*pow(cos(degree),2) - 30030*pow(cos(degree),4) + 90090*pow(cos(degree),6) -109395*pow((cos(degree)),8)+46189*pow(cos(degree),10))		
//Y011
					+pop[i][11]*(1/512.0)*sqrt(23/Pi)*(-693*pow(cos(degree),1) +15015*pow(cos(degree),3) - 90090*pow(cos(degree),5) +218790*pow((cos(degree)),7)-230945*pow(cos(degree),9)+88179*pow(cos(degree),11))
//Y012
					+pop[i][12]*(1/2048.0)*sqrt(25/Pi)*(231 -18018*pow(cos(degree),2) +225225*pow(cos(degree),4) - 1021020*pow(cos(degree),6) +2078505*pow((cos(degree)),8)-1939938*pow(cos(degree),10)+676039*pow(cos(degree),12))
					;
			}
			
			testScores[i] = Sum;
		}
/*		
		for (int k = 0; k <= PopMAX - 1; k++){
			cout << "Initial Score for " << k << ":" << testScores[k] << endl;
		}
*/
// The following matrix will hold the next generation's species. In the end of the evolution, we will make pop = nextPop, so the next evolution afterwards will act on nextPop 
		double nextPop[PopMAX][SphHarMAX] = {};
		
// Evolution Algorithim 1: Take the 10 species with the best scores and pass them onto nextPop
		double bestScores[10];
		int bestScoreIndicies[10];
		
		for (int i = 0; i <= 9; i++){
			bestScores[i] = testScores[i];
			bestScoreIndicies[i] = i;
		}
		
		for (int i = 10; i <= PopMAX - 1; i++){
			
			double minScore = bestScores[0];
			int minScoreIndex = 0;
			for (int j = 1; j <= 9; j++){
				if (bestScores[j] < minScore){
					minScore = bestScores[j];
					minScoreIndex = j;
				}
			}
			
			if (testScores[i] > minScore){
				bestScores[minScoreIndex] = testScores[i];
				bestScoreIndicies[minScoreIndex] = i;
			}
		}
		
		for (int i = 0; i <= 9; i++){
			cout << bestScoreIndicies[i] << ":" << bestScores[i] << endl;
		}
		cout << nextPop[10][0];
		for (int i = 0; i <= 9; i++){
			for (int j = 0; j <= SphHarMAX; j++){
				nextPop[i][j] = pop[bestScoreIndicies[i]][j];
			}
		}
		
		cout << "Generation after first Algorithim:" << endl;
		for (int i = 0; i <= 99; i++){
			for (int j = 0; j <= 12; j++){
				cout << nextPop[i][j] << "  ";
			}
			cout << endl;
		}
		
		cout << nextPop[10][0];
	}



/*double thetaM; // M for mean value
double phiM;
double stdev;
const double PI = 3.14159;
double norm = 0;
double gain[36][72];

// Enter the direction and standard deviation of the signal
cout << "Enter the direction and standard deviation (all in degrees) in the following format: theta phi stdev" << endl;
cin >> thetaM;
cin >> phiM,
cin >> stdev;
//double gaindB[180][360];

// We take each point on the sphere, find out its distance from (thetaM, phiM), and assign it a gain value based on the value of a normal PDF at that distance

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		double dist = distSph(i*5, j*5, thetaM, phiM) / 2;
		gain[i][j] = normPdf(dist, 0, stdev/ 72);
		//gaindB[i][j] = 10 * log10( normPdf(dist, 0, stdev/360) / ( norm /(360 * 180) ));
	}
}

// Now we normalize it in two steps, first calculate the sum over all points and then divide all value so that the sum  over all points equals 2592 (2592 = 36 * 72)

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		norm += gain[i][j] * sin(5* i * PI / 180.0)* 2 * PI; // We mulitply by the circumference at each point to normalize correctly
	}
}
cout << norm << endl;

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		gain[i][j] = 2592 * gain[i][j] / norm;
	}
}

// Finally, we output these values in a Mathematica-Friendly format 

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		cout << setprecision(10) << fixed <<"{" << i * 5 * PI / 180.0  << ", " << j * 5 * PI / 180.0 << ", " << gain[i][j] << "}, ";
	}
}

int end = 0;
cout << "Exit? (enter 1)" << endl;
cin >> end;*/

return 0;
}
