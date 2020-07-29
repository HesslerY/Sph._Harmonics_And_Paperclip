
// Written by: Suren Gourapura
// Date: 10/28/17
// Goal:     The goal is to evolve paperclip antennas using rotations as the genetics
//             These antenna will be maximizing curlyness about the z axis

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;

double Pi = 3.14159265359;
const int PopMAX = 100;
default_random_engine generator(time(NULL)); // add a seed if something is going wrong to watch for consistent results

// Code modified from https://codereview.stackexchange.com/questions/110793/insertion-sort-in-c
void insertionSort(double array[], int length){         // Sort array into greatest -> least
    int i,j;
    for (i = 1; i < length; i++) {
        double temp = array[i];
        for (j = i; j > 0 && array[j - 1] < temp; j--) {
            array[j] = array[j - 1];
        }
        array[j] = temp;
    }
}

void CoordTransform(double oldVec[], double rotx, double roty, double rotz, double newVec[]){
    // We are calculating the next unit vector using the previous unit vector and rotations
    double x = oldVec[0];
    double y = oldVec[1];
    double z = oldVec[2];
    // The below was calculated on mathematica using the 3D rotation matrices found here: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
    // We are performing rotx first, then roty, then rotz
    newVec[0] = sin(rotx)*(y*sin(roty)*cos(rotz) + z*sin(rotz)) + cos(rotx)*(z*sin(roty)*cos(rotz) - y*sin(rotz)) + x*cos(roty)*cos(rotz);
    newVec[1] = sin(rotz)*(y*sin(rotx)*sin(roty) + x*cos(roty)) + cos(rotx)*(z*sin(roty)*sin(rotz) + y*cos(rotz)) - z*sin(rotx)*cos(rotz);
    newVec[2] = y*sin(rotx)*cos(roty) + z*cos(rotx)*cos(roty) - x*sin(roty);
}

void CrossProduct(double oldVec[], double newVec[], double crossVec[]){
    double x = oldVec[0];
    double y = oldVec[1];
    double z = oldVec[2];
    
    double a = newVec[0];
    double b = newVec[1];
    double c = newVec[2];
    
    crossVec[0] = c*y - b*z;
    crossVec[1] = a*z - c*x;
    crossVec[2] = b*x - a*y;
}

double FScore(int numSeg, double rotx[], double roty[], double rotz[]){
    // First, we need to convert from rotations to unit vector coordinates.
    // Each unit vector corresponds to the direction that line segment is pointing relative to fixed coordinates
    // Initialize the converted array
    double unitVecs[numSeg][3];
    
    // The first node is rotated from vertical (0,0,1)
    double newVec[3];
    double oldVec[] = {0,0,1};
    CoordTransform(oldVec, rotx[0], roty[0], rotz[0], newVec);
    for (int i = 0; i < 3; i++){
        unitVecs[0][i] = newVec[i];
    }
    
    // Now, we use loops to generate the rest of the coordinates
    for (int i = 1; i < numSeg; i++){
        for (int j = 0; j < 3; j++){                                    // Get the old vector (the point before point i)
            oldVec[j] = unitVecs[i - 1][j];
        }
        CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);    // Calculate the new vector based on old vector and rotations
        for (int j = 0; j < 3; j++){
            unitVecs[i][j] = newVec[j];                                // Add new vectors onto cartesian
        }
    }
    
    // With [numSeg] unit vectors, we can calculate [numSeg-1] cross product vectors.
    // The magnitude of these vectors in the z direction (arbitrary choice) gives us the fitness score.
    double crossVec[numSeg - 1][3];
    
    for (int i = 0; i < numSeg - 1; i++){
        for (int j = 0; j < 3; j++){
            oldVec[j] = unitVecs[i][j];        // Initialize the old vector
        }
        for (int j = 0; j < 3; j++){
            newVec[j] = unitVecs[i + 1][j];    // Initialize the new vector
        }
        CrossProduct(oldVec, newVec, crossVec[i]);
    }
    
    // Simply sum the z component of the cross product vectors to get the fitness score
    double fScore = 0;
    
    for (int i = 0; i < numSeg - 1; i++){
        fScore += crossVec[i][2];
    }

    if(fScore < 0.0)
      {
	fScore = fScore*(-1.0);
      }
    
    return fScore;
}

double RotToCartesian(int numSeg, double rotx[], double roty[], double rotz[], double xcoord[], double ycoord[], double zcoord[]){
    // First, we need to convert from rotations to unit vector coordinates.
    // Each unit vector corresponds to the direction that line segment is pointing relative to fixed coordinates
    // Initialize the converted array
    double unitVecs[numSeg][3];
    
    // The first node is rotated from vertical (0,0,1)
    double newVec[3];
    double oldVec[] = {0,0,1};
    CoordTransform(oldVec, rotx[0], roty[0], rotz[0], newVec);
    for (int i = 0; i < 3; i++){
        unitVecs[0][i] = newVec[i];
    }
    
    // Now, we use loops to generate the rest of the coordinates
    for (int i = 1; i < numSeg; i++){
        for (int j = 0; j < 3; j++){                                    // Get the old vector (the point before point i)
            oldVec[j] = unitVecs[i - 1][j];
        }
        CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);    // Calculate the new vector based on old vector and rotations
        for (int j = 0; j < 3; j++){
            unitVecs[i][j] = newVec[j];                                // Add new vectors onto cartesian
        }
    }
    
    // Now, we need to convert the unit vectors into the actual coordinates. The first point is {0,0,0} and the second point is the same as the first unit vector
    // Additional coordinates are made by adding the unit vector onto the previous coordinate
    xcoord[0] = 0;
    ycoord[0] = 0;
    zcoord[0] = 0;
    
    for (int i = 0; i < numSeg + 1; i++){ // I honestly have no idea why it is numSeg + 1 instead of numSeg, but the latter doesn't work!
        xcoord[i+1] = xcoord[i] + unitVecs[i][0];
        ycoord[i+1] = ycoord[i] + unitVecs[i][1];
        zcoord[i+1] = zcoord[i] + unitVecs[i][2];
    }
    return 0;
}

// Ryan's test functions
int roulette(double fitness[], int Pop); // roulette selection function

int tournament(double fitness[], int Pop); // tournament selection function

int rank(double fitness[], int Pop); // rank tournament selection function

vector<vector<vector<double> > > crossover(vector<vector<vector<double> > > & cross);

vector<vector<vector<double> > > simple_mutation(vector<vector<vector<double> > > & mutation, double mut_chance, double std_dev);

void Results(double fitness[], int generations, double mut_chance, string run_num, vector<double> high_score, vector<double> average_score, int Tour, int Roul);

void data(int generations, string run_num, vector<double> high_score, vector<double> average_score);
//

int main()
{
    const int numSeg = 10;
    int Gen = 0;
    srand(time(NULL));
    
    // Enter the number of antenna segments
    cout << "Enter the number of generations:" << endl;
    // cin >> numSeg;
    cin >> Gen;

    //Ryan's on or off switches for tournament and roulette
    int Roul;
    int roul_num;
    int tour_num;
    int Tour;
    double mut_chance;
    double std_dev;
    string run_num;
    cout << "Enter run number: ";
    cin >> run_num;
    cout << "Would you like to use Ryan's Algorithm? (1 = yes)  " ;
    cin >> Roul;
    if (Roul == 1)
      {
	Tour = Roul;
      }
    if(Roul == 1 && Tour == 1)
      {
	cout << "Please enter the desired mutation probability in decimal notation (between 0.0-1.0): " ;
	cin >> mut_chance;
	while(mut_chance > 1.0 || mut_chance < 0.0)
	  {
	    cout << "Enter a valid mutation probability: ";
	    cin >> mut_chance;
	  }
	cout << "Please enter the (positive) standard deviation: ";
	cin >> std_dev;
      }
    if(Roul == 1 && Tour == 1)
      {
	cout << "Please enter how many individuals will be sent to roulette (10 max). The rest will be sent to tournament: ";
	cin >> roul_num;
	while( roul_num < 0 || roul_num > 10)
	  {
	    cout << "Enter a valid integer: ";
	    cin >> roul_num;
	  }
	tour_num = 10 - roul_num;
      }
    
    
    // Create the population with user specified number of line segments
    double pop[PopMAX][numSeg][3];
    double nextPop[PopMAX][numSeg][3];
    
                                                                                                                               
    // Randomly fill the population. First step is to give each x, y, and z rotation a random value between 0 and 2 pi.
    for (int i = 0; i < PopMAX; i++){                                       // for each antenna
        for (int j = 0; j < numSeg; j++){                                      // for each segment
            for (int k = 0; k < 3; k++){                                      // for each x, y, and z rotation
                pop[i][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi;     // gives a random value between 0 and 2 pi
            }
        }
    }


    // create vectors for storing generation information
    vector<double> high_score;
    vector<double> average_score;

    
    
    // Now, we evolve. We loop the following code for each generation
    for (int g = 1; g <= Gen; g++){
        //    cout << endl << "Generation: " << g << endl;
        
        // The first objective is to calculate how well our current population's species are doing
        double testScores[PopMAX] = {};
        
        // To test the curent population, we first reorganize the 2D array (numSeg x 3) into three 1D arrays, each containing all values for one rotation
        double rotx[numSeg], roty[numSeg], rotz[numSeg];
        
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < numSeg; j++){
                rotx[j] = pop[i][j][0];
                roty[j] = pop[i][j][1];
                rotz[j] = pop[i][j][2];
            }
            testScores[i] = FScore(numSeg, rotx, roty, rotz); // Each antenna gets its score recorded. The rot arrays are rewritten for each species
        }
        
        // The following matrix will hold the next generation's species. In the end of the evolution, we will make pop = nextPop, so the next evolution afterwards will act on nextPop
	// double nextPop[PopMAX][numSeg][3] = {};
        
	// store information at the begining of each generation
	double rank_scores[PopMAX] = {};
	double total_score=0;
	for(int x=0; x<PopMAX; x++)
	  {
	    rank_scores[x] = testScores[x];
	    total_score = total_score + testScores[x];
	  }

	insertionSort(rank_scores, PopMAX);
	high_score.push_back(rank_scores[0]);
	average_score.push_back((total_score/PopMAX));



	if(Roul != 1 && Tour != 1)
	  {
        
        // Evolution Algorithim 1: Take the 10 species with the best scores and pass them onto nextPop
        // reproduction
        // We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop
        // Also, a temporary ordered testScores array: rankedTestScores
        double rankedPop[PopMAX][numSeg][3] = {};
        double rankedTestScores[PopMAX] = {};
        
        // First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
        for (int i = 0; i < PopMAX; i++){
            rankedTestScores[i] = testScores[i];
        }
        
        insertionSort(rankedTestScores, PopMAX);
        
        // Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < PopMAX; j++){
                if (testScores[j] == rankedTestScores[i]){
                    
                    for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
                        for (int l = 0; l < 3; l++){
                            rankedPop[i][k][l] = pop[j][k][l];
                        }
                    }
                    
                }
            }
        }
        
        // We print out the highest ranking species's score and it's array in mathematica format, for easy plotting
        cout << rankedTestScores[0] << endl;
        
        // Finally, we copy over the top 10 best species from pop to nextPop
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < numSeg; j++){
                for (int k = 0; k < 3; k++){
                    nextPop[i][j][k] = rankedPop[i][j][k];
                }
            }
        }
        
        
        
        
        cout << "Generation " << g << " top 5:" << endl;
        for (int i = 0; i < 5; i++){
            cout << "# " << i+1 << " With score "<< rankedTestScores[i]<< " : " << endl;
            for (int j = 0; j < numSeg; j++){
                cout << "Node " << j << " : X-Rot = " << rankedPop[i][j][0];
                cout << ", Y-Rot = " << rankedPop[i][j][1];
                cout << ", Z-Rot = "<< rankedPop[i][j][2] << endl;
            }
            cout << endl;
        }
        cout << "# " << 100 << " With score "<< rankedTestScores[99]<< " : " << endl;
        for (int j = 0; j < numSeg; j++){
            cout << "Node " << j << " : X-Rot = " << rankedPop[99][j][0];
            cout << ", Y-Rot = " << rankedPop[99][j][1];
            cout << ", Z-Rot = "<< rankedPop[99][j][2] << endl;
        }
	  
  
        // Evolution Algorithim 2: Take 10 random species, find the one with the best score, and [randomly mutate one of it's rotations to obtain an offspring] 10 times.
        // Do this whole algorithim 3 times
	// mutation
        for(int a2 = 1; a2 <= 4; a2++){
            
            int choose10[10]; // Create an array with 10 random values, 0-99. This array determines the 10 random species that will undergo a tournament selection (a.k.a. simply choosing the highest score species out of the 10)
            for (int i = 0; i < 10; i++){
                choose10[i]= rand() % 100;
            }
            
            // Since we have pop already organized from best to worst, we simply find the lowest value in choose10 to find the winner of the tournament. The winning species is called: Algorithim 2 Best Value
            int Alg2BestVal = choose10[0]; // Assume the best species is the first one
            
            for (int i = 1; i < 10; i++){
                if (Alg2BestVal > choose10[i]){  // If another species in the choose10 array has a lower value, it is now the best species
                    Alg2BestVal = choose10[i];
                }
            }
            
            // Now, we mutate one part of the species's 2D array, pop[Alg2BestVal][][], to create an offspring. We do this 10 times.
            int mutateLocation[10][2]; // Create the locations for mutation of the 10 offspring. Each location is 2D: which node, which rotation
            for (int i = 0; i < 10; i++){
                mutateLocation[i][0] = rand() % numSeg;    // What segment will be mutated
                mutateLocation[i][1] = rand() % 3;        // What rotation will be mutated
            }
            
            // We do the process below 10 times, (to spots 10*a2 to 19*a2 in nextPop)
            for (int i = 0; i < 10; i++){
                for (int j = 0; j < numSeg; j++){
                    for (int k = 0; k < 3; k++){
                        nextPop[i + 10*a2][j][k] = rankedPop[Alg2BestVal][j][k];    // Copy the Alg2BestVal species over. Note: a2 is added to be able to do the whole of Algorithm 2, 3 times
                    }
                }
                for (int j = 0; j < numSeg; j++){
                    if (mutateLocation[i][0] == j){            // If the node location is right
                        for (int k = 0; k < 3; k++){
                            if (mutateLocation[i][1]==k){    // If the rotation location is right
                                nextPop[i + 10*a2][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi;  // Mutate this location
                            }
                        }
                    }
                }
            }
        }
        
        
        
        // Evolution Algorithim 3: Take 20 random species, run two seperate tournaments to find 2 parents, and [swap a random array location with each other to obtain two offspring (for both combinations)] 5 times.
        // We do this Algorithm 5 times, for a total of 50 offspring
	//crossover
        for(int a3 = 1; a3 <= 5; a3++){
            
            int choose10A[10], choose10B[10]; // Create 2 arrays with 10 random values, 0-99. These arrays determine the 20 random species that will undergo two seperate tournament selections
            int Alg3BestValA = 0, Alg3BestValB = 0;
            
            while (Alg3BestValA == Alg3BestValB){ // To make sure that the 2 parents aren't the same species, we run the following as long as they are the same
                for (int i = 0; i < 10; i++){
                    choose10A[i]= rand() % 100;
                    choose10B[i]= rand() % 100;
                }
                
                // Since we have pop already organized from best to worst, we simply find the lowest value in choose10 to find the winners of the tournaments. The winning species are called: Algorithim 2 Best Value A/B
                Alg3BestValA = choose10A[0]; // Assume the best species are the first ones
                Alg3BestValB = choose10B[0];
                for (int i = 1; i < 10; i++){
                    if (Alg3BestValA > choose10A[i]){  // If a lower value is found, make best value that lower value
                        Alg3BestValA = choose10A[i];
                    }
                    if (Alg3BestValB > choose10B[i]){
                        Alg3BestValB = choose10B[i];
                    }
                }
            }
            
            // Now, we swap one part of the species's arrays in both ways to create 2 offspring. We do this 5 times, normalizing after each one
            int swapLocation[5]; // Create the locations for swapping
            for (int i = 0; i < 5; i++){
                swapLocation[i] = rand() %  100;
            }
            
            // We do the process below 5 times, (to spots 50*a3 to 59*a3 in nextPop)
            for (int i = 0; i < 5; i++){
                // First, we copy over the two parents in the 10 offspring spots in nextPop, in an A,B,A,B,... pattern
                for (int j = 0; j < numSeg ; j++){
                    for (int k = 0; k<3; k++){
                        nextPop[(i*2) + 10*a3 + 40][j][k]= rankedPop[Alg3BestValA][j][k];    // Copy the Alg3BestVal[A and B] species over. Note: a3 is added to be able to do the whole of Algorithm 3, 5 times
                        nextPop[(i*2 + 1) + 10*a3 + 40][j][k]= rankedPop[Alg3BestValB][j][k];
                    }
                }
                for (int j = 0; j < numSeg ; j++){
                    for (int k = 0; k<3; k++){
                        if(swapLocation[i] == j){  // make sure that the location to be swapped is chosen using swapLocation[]
                            double temp = nextPop[(i*2) + 10*a3 + 30][j][k];
                            nextPop[(i*2) + 10*a3 + 40][j][k] = nextPop[(i*2+1) + 10*a3 + 30][j][k];
                            nextPop[(i*2+1) + 10*a3 + 40][j][k] = temp;
                        }
                    }
                }
                
            }
        }
      
	  }


       


	// Ryan's algorithm
	if(Tour == 1 && Roul == 1)
	  { 
	    
	    //pre-select individuals and set them into arrays for later use
	    vector<int> tour_indiv;                                                                                                                                                                                                                     vector<int> roul_indiv;     
	    for(int i=0; i<PopMAX; i++)
	      {
		if(Tour == 1)
		  {
		    tour_indiv.push_back(tournament(testScores, PopMAX));
		  }
		if(Roul ==1)
		  {
		    roul_indiv.push_back(roulette(testScores, PopMAX));
		  }
	      }
	  
	    
	    // Reproduction
	    
	    for(int i=0; i<10; i++)
	      {
		for(int j=0; j<numSeg; j++)
		  {
		    for(int k=0; k<3; k++)
		      {
			if(i < roul_num)
			  {
			    nextPop[i][j][k] = pop[roul_indiv[i]][j][k];
			  }
			else
			  {
			    nextPop[i][j][k] = pop[tour_indiv[i]][j][k];
			  }
		      }
		  }
	      }
	    
	    // crossover
	    // prime vectors to store information
	    vector<vector<vector<double> > > cross;
	    vector<vector<vector<double> > > children;
	    
	    for(int i=0; i<90; i++)
	      {
		for(int j=0; j<numSeg; j++)
		  {
		    for(int k=0; k<3; k++)
		      {
			if (i<2)
			  {
			    cross.push_back(vector<vector<double> > ());
			    cross[i].push_back(vector<double> ());
			    cross[i][j].push_back(0);
			  }
			children.push_back(vector<vector<double> > ());
			children[i].push_back(vector<double> ());
			children[i][j].push_back(0);
		      }
		  }
	      }
  
	    for(int i=0; i<90; i=i+2)
	      {
		for(int j=0; j<10; j++)
		  {
		    for(int k=0; k<3; k++)
		      {
			if ( i<(9*roul_num))
			  {
			    cross[0][j][k] = pop[roul_indiv[i+roul_num]][j][k];
			    cross[1][j][k] = pop[roul_indiv[i+1+roul_num]][j][k];
			  }
			else
			  {
			    cross[0][j][k] = pop[tour_indiv[i+tour_num]][j][k];
			    cross[1][j][k] = pop[tour_indiv[i+1+tour_num]][j][k];
			  }
		      }
		  }
		cross = crossover(cross);

		for(int x=0; x<10; x++)
		  {
		    for( int y=0; y<3; y++)
		      {
			children[i][x][y] = cross[0][x][y];
			children[i+1][x][y] = cross[1][x][y];
		      }
		  }
	      }
		 

	    
	    // mutations
	    
	    children = simple_mutation(children, mut_chance, std_dev);

	    for(int x=10; x<100; x++)
	      {
		for(int y=0; y<10; y++)
		  {
		    for( int z=0; z<3; z++)
		      {
			nextPop[x][y][z] = children[x-10][y][z];
		      }
		  }
	      }
	  }
	




 
        /*
         // Evolution Algorithim 4: Introduce 10 random species into the population
         // First step is to give each value in each species random number between -1 and 1. Then, we normalize the species values to 1
         for (int i = 0; i < 10; i++){
         for (int j = 0; j < SphHarMAX; j++){
         nextPop[90 + i][j] = ((double)rand() / (double)(RAND_MAX))*2- 1; // gives a random value between -1 and 1
         }
         }
         
         // Now we normalize nextPop: We divide all coefficients by their integral over all space to normalize to 1
         integral = 0;
         for(int i =0; i < PopMAX; i++){
         for(int j =0; j < SphHarMAX; j++){
         integral += pow(nextPop[i][j], 2);
         }
         for(int j =0; j < SphHarMAX; j++){
         nextPop[i][j] = nextPop[i][j] / sqrt(integral);
         }
         integ         }
         */
        // Finally, we equate the old pop to the new pop, allowing the next loop to operate on the new population
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < numSeg; j++){
                for (int k = 0; k < 3; k++){
                    pop[i][j][k] = nextPop[i][j][k];
                }
            }
        }
        
    }    // End of Evolution
    

    
    // The work is done, now time to see the results of the final generation! We follow the same ranking protocol
    double finalScores[100];
    
    
    // Now, to test the curent population
    double rotx[numSeg], roty[numSeg], rotz[numSeg];
    
    for (int i = 0; i < PopMAX; i++){
        for (int j = 0; j < numSeg; j++){
            rotx[j] = pop[i][j][0];
            roty[j] = pop[i][j][1];
            rotz[j] = pop[i][j][2];
        }
        finalScores[i] = FScore(numSeg, rotx, roty, rotz); // Each antenna gets its score recorded. The rot arrays are rewritten for each species
    }
    
    // We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop
    // Also, a temporary ordered testScores array: rankedTestScores
    double rankedPop[PopMAX][numSeg][3] = {};
    double rankedFinalScores[PopMAX] = {};
    // First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
    for (int i=0; i < PopMAX; i++){
        rankedFinalScores[i] = finalScores[i];
    }
    
    insertionSort(rankedFinalScores, PopMAX);
    
    // Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
    for (int i = 0; i < PopMAX; i++){
        for (int j = 0; j < PopMAX; j++){
            if (finalScores[j] == rankedFinalScores[i]){
                
                for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
                    for (int l = 0; l < 3; l++){
                        rankedPop[i][k][l] = pop[j][k][l];
                    }
                }
                
            }
        }
    }
    /*
     for (int i = 0; i < PopMAX; i++){
     for (int k = 0; k < PopMAX; k++){
     if (finalScores[k] == rankedFinalScores[i]){
     for (int l = 0; l < SphHarMAX; l++){ // We need to copy the whole of the species in the kth position to rankedPop in the ith position
     rankedPop[i][l] = pop[k][l];
     }
     
     }
     }
     } */
    
    // We print out the highest ranking species's scores and it's arrays in mathematica format, for easy plotting
    cout << endl << endl << "Final Results:" << endl;
    
    for (int i = 0; i < 5; i++){
        cout << "# " << i+1 << " With score "<< rankedFinalScores[i]<< " : " << endl;
        for (int j = 0; j < numSeg; j++){
            cout << "Node " << j << " : X-Rot = " << rankedPop[i][j][0];
            cout << ", Y-Rot = " << rankedPop[i][j][1];
            cout << ", Z-Rot = "<< rankedPop[i][j][2] << endl;
        }
        cout << endl;
    }
    
    // To print in Mathematica friendly format, we first need to convert the rotations into unit vectors, then into cartesian coordinates
    
    
    
    cout << endl << endl << "Final Results in Mathematica Format:" << endl;
    
    for (int i = 0; i < 5; i++){
        cout << "# " << i+1 << " With score "<< rankedFinalScores[i]<< " : " << endl;
        double rotx[numSeg], roty[numSeg], rotz[numSeg];
        double xcoord[numSeg+1], ycoord[numSeg+1], zcoord[numSeg+1];
        
        for (int j = 0; j < numSeg; j++){
            rotx[j] = rankedPop[i][j][0];
            roty[j] = rankedPop[i][j][1];
            rotz[j] = rankedPop[i][j][2];
        }
        
        RotToCartesian(numSeg, rotx, roty, rotz, xcoord, ycoord, zcoord);
        
        cout << "line = Line [{";
        for (int i = 0; i < numSeg; i++){
            cout << "{" << xcoord[i] << ", " << ycoord[i] << ", "<< zcoord[i] << "}, ";
        }
        cout << "{" << xcoord[numSeg] << ", " << ycoord[numSeg] << ", "<< zcoord[numSeg] << "}}] " << endl;
    }

    Results(rankedFinalScores, Gen, mut_chance, run_num, high_score, average_score, Tour, Roul);
    data(Gen, run_num, high_score, average_score);
  
    return 0;
}


int roulette(double fitness[], int Pop)
{
  // re-write fitness scores to ignore negative fitness scores and find the sum of positive fitness scores
  vector<double> adjusted_fit;
  double total;
  double min =0;
  
  for(int x=0; x<100; x++)
    {
      if(fitness[x] < min)
	{
	  min = fitness[x];
	}
    }
	    
    
  for(int i=0; i<Pop; i++)
    {
      /*
      if (fitness[i] <= 0.0)
	{
	  i++;
	}/
	else */
	
	  adjusted_fit.push_back(fitness[i]);
	
      total = total + adjusted_fit[i];
    }
  
  
  // run a roulette selection
  int chosen;
  int check;
  double select;
  double over;
  double under = 0;
  uniform_real_distribution<double> choice(0.0, total);

  select = choice(generator);

  for (int x=0; x<adjusted_fit.size(); x++)
    {
      under = under + adjusted_fit[x];
      if( x == (adjusted_fit.size()-1))
	{
	  check = x;
	  x = adjusted_fit.size();
	}
      else
	{
	  over  = under + adjusted_fit[x+1];
	}
      
      if (select <= under)
	{
	  check  = x;
	  x = adjusted_fit.size();
	}
    }

  for( int y=0; y<Pop; y++)
    {
      if(fitness[y] == adjusted_fit[check])
	{
	  chosen = y;
	}
    }
  // return the selected individual     

  return(chosen);
}


int tournament(double fitness[], int Pop)
{
  int chosen;
  double select;
  vector<double> tourney;
  uniform_real_distribution<double> choice(0.0, 99.0);
  
  // randomly select 10 individuals from the population
  for (int i=0; i<10; i++)
    {
      select = choice(generator);
      for (int j=0; j<100; j++)
	{ 
	  if ( j == 99)
	    {
	      tourney.push_back(fitness[j]);
	      i = i+1;
	    }
	  else if(select == j || select < (j+1) && select > j)
	    {
	      tourney.push_back(fitness[j]);
	      i=i+1;
	    }
	}
    }
  double max = tourney[0];
  
  // find the highest fitness score in that group
  for (int k=0; k<10; k++)
    {
      if (tourney[k] > max)
	{
	  max = tourney[k];
	}
    }

  // find the individual the score belongs to and return the individual
  for(int x=0; x<100; x++)
    {
      if( fitness[x] == max)
	{
	  chosen = x;
	  x = 100;
	}
      if( x == 99)
	{
	  chosen = 99;
	  x = 100;
	}
    }
  return(chosen);
}

int rank(double fitness[], int Pop) // this is just for Ryan to have in case it wants to be used in the future
{
  int chosen;
  return (chosen);
}



vector<vector<vector<double> > > crossover(vector<vector<vector<double> > > & cross)
{
  double swap;
  uniform_real_distribution<double> choice(0.0, 1.0);
  for(int i=0; i<10; i++)
    {
      for (int j=0; j<3; j++)
	{
	  swap = choice(generator);
	  if(swap < 0.5)
	    {
	      double temp = cross[0][i][j];
	      cross[0][i][j] = cross[1][i][j];
	      cross[1][i][j] = temp;
	    }
	}
    }
	  
  return(cross);
}

vector<vector<vector<double> > > simple_mutation(vector<vector<vector<double> > > & mutation, double mut_chance, double std_dev)
{
  
  double mut;
  double prob;
  //default_random_engine generator(1);
  uniform_real_distribution<double> chance(0.0, 1.0);
  // normal_distribution<double> angle(0.0, 2*M_PI);

  for(int i=0; i<90; i++)
    {
      for(int j=0; j<10; j++)
	{
	  for(int k=0; k<3; k++)
	    {
	      prob = chance(generator);
	    
	      if(prob <= mut_chance)
		{
		  double mean = mutation[i][j][k];
		  normal_distribution<double> angle(mean, std_dev);
		  
		  mut = angle(generator);
		  while(mut <0 || mut > 2*M_PI)
		    {
		      if(mut < 0) 
			{
			  mut = mut + 2.0*M_PI; 
			}
		      else if(mut > 2*M_PI)
			{
			  mut= mut - 2*M_PI;
			}
		    }
		  mutation[i][j][k] = mut;
		}
	    }
	}
    }
  
  return(mutation);
		  
}

void Results(double fitness[], int generations, double mut_chance, string run_num, vector<double> high_score, vector<double> average_score, int Tour, int Roul)
{
  double total=0;
  double average;
  double maximum=0 ;
  int max;

  for(int m=0; m<100; m++)
    {
      if (high_score[m] > maximum)
	{
	  maximum= high_score[m];
	  max = m;
	}
    }

  if(fitness[0] > maximum)
    {
      maximum = fitness[0];
      max = generations +1;
    }
  for (int i=0; i<100; i++)
    {
      total = total + fitness[i];
    }
  average = total/100.0;
  ofstream Run;
  Run.open("Run_"+ run_num + ".txt");
  Run <<"This is a test of paperclips1.0.1.cpp \n";
  if(Tour == 1)
    {
      Run << "USING TOURNAMENT: "<< endl;
    }
  if(Roul == 1)
    {
      Run << "USING ROULETTE: " << endl;
    }
  Run <<"GENERATIONS: "<< generations << endl;
  Run <<"MUTATION PROBABILITY: " << mut_chance << endl;
  Run <<"HIGHEST FITNESS SCORE: " << maximum << " IN GEN: " << max << endl;
  Run <<"TOTAL COMBINED FITNESS SCORE: "<< total << endl;
  Run <<"AVERAGE FITNESS SCORE: " << average << endl;

  Run <<"GENERATION TRENDS: " << endl;
  
  for (int x=0; x<generations; x++)
    {
      Run << "Gen " << x << ": " << endl;
      Run << "      High Score: " << high_score[x] << endl;
      Run << "      Ave. Score: " << average_score[x] << endl;
      Run << endl;
    }

  Run << endl << "Final Generation Individuals Fitness Scores: "<< endl;

  for(int j=0; j<100; j++)
    {
      Run <<"INDIVIDUAL "<< j+1 << " : " << fitness[j] << endl;
    }

  Run.close();
}

void data(int generations, string run_num, vector<double> high_score, vector<double> average_score)
{
  ofstream Run;
  
  Run.open("Run_"+run_num+".csv");

  Run << "This is a csv file for generations. See the same numbered Run.txt for other information." << endl;
  
  Run << "Generation: Gen High: Gen average:" << endl;

  for(int i=0; i<generations; i++)
    {
      Run << i << "," << high_score[i] << "," << average_score[i] << endl;
    }
  Run.close();
}
