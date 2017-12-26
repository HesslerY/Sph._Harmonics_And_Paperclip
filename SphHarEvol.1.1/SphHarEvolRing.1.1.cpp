
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
	
// Randomly fill the population. First step is to give each value in each species random number between -1/2 and 1/2. Then, we normalize the species values to 1
	for (int i = 0; i <= PopMAX - 1; i++){
		for (int j = 0; j <= SphHarMAX - 1; j++){
			pop[i][j] = ((double)rand() / (double)(RAND_MAX))-.5; // gives a random value between -1/2 and 1/2
		}
	}
	
// We add up the values for each species and divide them by the sum to normalize to 1
	for (int i = 0; i <= PopMAX - 1; i++){
		double counter = 0;
		for (int j = 0; j <= SphHarMAX - 1; j++){
			counter += fabs(pop[i][j]);
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
			double Min;
			double Max;
			
// Make sure that the maximum is less than or equal to 180 degrees and the minimum is greater than or equal to 0
			if (Theta + Spread/2.0 > 180){Max = 180;}
			else{Max = Theta + Spread/2.0;}
			
			if (Theta - Spread/2.0 < 0){Min = 0;}
			else{Min = Theta - Spread/2.0;}
			
// Sum the spherical harmonic values by degree increments from the Min to Max theta values
// Note: I would have loved to make this a seperate function, but I found that passing matrices is surprisingly challenging in c++. Sorry!
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
// Finally, we give the summed values of the spherical harmonic from theta min to theta max as the test score (a rough kind of integration by 1 degree increments)			
			testScores[i] = Sum;
		}

// The following matrix will hold the next generation's species. In the end of the evolution, we will make pop = nextPop, so the next evolution afterwards will act on nextPop 
		double nextPop[PopMAX][SphHarMAX] = {};
		
// Evolution Algorithim 1: Take the 10 species with the best scores and pass them onto nextPop

// We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop 
// Also, a temporary ordered testScores array: rankedTestScores
		double rankedPop[PopMAX][SphHarMAX] = {};
		double rankedTestScores[PopMAX] = {};
// First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
		for (int i=0; i <= PopMAX-1; i++){
			rankedTestScores[i] = testScores[i];
		}

		insertionSort(rankedTestScores, PopMAX);
		
// Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
		for (int i=0; i <=PopMAX-1; i++){
			for (int k=0; k <=PopMAX-1; k++){
				if (testScores[k] == rankedTestScores[i]){
					
					for (int l=0; l <=SphHarMAX-1; l++){ // We need to copy the whole of the species in the kth position to rankedPop in the ith position
						rankedPop[i][l] = pop[k][l];
					}
					
				}
			}
		} 
				
		cout << "First Score in Mathematica Format with score "<< rankedTestScores[0]<< " : " << endl;
		for (int i = 0; i <= SphHarMAX-1; i++){
			cout << "m1Y0" << i<< " = " << rankedPop[0][i] << endl;
		}
		
		for (int i = 0; i <= 10-1; i++){
			for (int j = 0; j <= SphHarMAX-1; j++){
				nextPop[i][j]= rankedPop[i][j];
			}

		}
		/*
		cout << "Ordered Generation after first Algorithim:" << endl;
		for (int i = 0; i <= PopMAX-1; i++){
			cout << "Test Score: " << rankedTestScores[i] << " ";
			for (int j = 0; j <= SphHarMAX-1; j++){
				cout << nextPop[i][j] << "  ";
			}
			cout << endl;
		}		
		*/
// Evolution Algorithim 2: Take 10 random species, find the one with the best score, and [randomly mutate one of it's array values to obtain an offspring]*10. Do this 3 times
		for(int a2 = 1; a2 <= 3; a2++){		
		
			int choose10[10];
			for (int i = 0; i <= 9; i++){
				choose10[i]= rand() % 100;
			}
		
// Since we have pop already organized from best to worst, we simply find the lowest value in choose10
			int Alg2BestVal = choose10[0];
		
			for (int i = 1; i <= 9; i++){
				if (Alg2BestVal > choose10[i]){
					Alg2BestVal = choose10[i];		
				}
			}
// Now, we mutate one part of pop[Alg2BestVal][] to create an offspring. We do this 10 times, normalizing after each one
			int mutateLocation[10]; // Create the location for mutation of the 10 offspring
			for (int i = 0; i <= 9; i++){
				mutateLocation[i]= rand() % SphHarMAX;
			}
// We do the process below 10 times, (to spots 10-19 in nextPop)
			for (int i = 0; i <= 9; i++){
			
				for (int j = 0; j <= SphHarMAX-1; j++){
					nextPop[i + 10*a2][j]= rankedPop[Alg2BestVal][j];	// Copy the Alg2BestVal species over
				}
				for (int j = 0; j <= SphHarMAX-1; j++){
					if(mutateLocation[i] == j){  // make sure that the location to be mutated is chosen randomly, by using mutateLocation[]
						nextPop[i + 10*a2][j] = ((double)rand() / (double)(RAND_MAX))-.5;// Mutate this location	
					}
				}	
// Now, we normalize this mutated species
				double counter = 0;
				for (int j = 0; j <= SphHarMAX - 1; j++){
					counter += fabs(nextPop[i + 10*a2][j]);
				}
				for (int j = 0; j <= SphHarMAX - 1; j++){
					nextPop[i + 10][j] = nextPop[i + 10*a2][j] / counter;
				}				
			}		
			

					
		}
/*		
		cout << "Generation after Second Algorithim:" << endl;
		for (int i = 0; i <= 39; i++){ //FIX
			for (int j = 0; j <= SphHarMAX-1; j++){
				cout << nextPop[i][j] << "  ";
			}
			cout << endl;
		}
*/			
// Finally, we equate the old pop to the new pop, allowing the next loop to operate on the new population
		for (int i = 0; i <= PopMAX-1; i++){
			for (int j = 0; j <= SphHarMAX-1; j++){
			pop[i][j] = nextPop[i][j];
			}
		}			 
	}	
		
		
		
			






return 0;
}
