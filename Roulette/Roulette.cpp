// Welcome to the new Paperclips Roulette Algorithm//
// Our goal is to test this algorithm (based on the bicone's roulette algorithm) against the existing tournament
// algorithm so see what differences might arise and see if it could be usefull to the current bicone evolution

// Code Starts Here
// **IMPORTANT** compile using g++ -std=c++11 Roulette.cpp -o Roulette.exe

#include <time.h>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <thread>

using namespace std;

// Headers

void dataRead(vector<vector<vector<float>>> &varInput, vecotr<float> &fitness);

void dataWrite(int numChildren, vector<vector<vector<float>>> &varVector, int freq_coeffs, vector<double> freqVector);

void checkConvergence(vector<vector<vector<float>>> &varInput, vector<float> &fitness);

void roulette(vector<vector<vector<float>>> &varInput, vector<vector<vector<float>>> &varOutput, vecotr<float> &fitness);


//
// Global Constants

double Min_Freq = .08333;
double Max_Freq = 1.0667;
double Freq_Step = 0.01667;

//
// DNA Constants

const int NSections = 1;
const int NVars = 3;
const int Partent_NO = 2;
const int DNA_Garbage_End = 9;

//
// Algorithm Constants

int NPop;
const float Mutability =0.6f;

//
// Statistic constants

const float Convergence =0.00;

//
// Genes
// add to this once we find all of the genes for the paperclips 



// Main Function

int main(int argc, char const *argv[])
{

  double Geoscale_Factor = stod(argv[3]);

  // genes things
  //
  //
  //
  //

  // frequencies scaled inversely with dimensions

  Max_Freq *= Geoscale_Factor;
  
  Min_Freq *= Geoscale_Factor;

  Freq_Step *= Geoscale_Factor;


  // Define NPop

  NPop = atoi(argv[2]);
 
  vector<vector<vector<float>>> varInput (NPop,vector<vector<float> >(Nsections, vector <float>(NVars, 0.0f)));

  vector<float> fitness (NPop, 0.0f);

  vector<vector<vector<float>>> varOutput (NPop, vector<vector<float> >(NSections, vector <float>(Nvars, 0,0f)));

  int freq_coeffs = round((Max_Freq - Min_Freq) / Freq_Step + 1);

  vector<double> freqVector (freq_coeffs, 0.0);

  freqVector[0] = Min_Freq;

  for(int i=1; i < freq_coeffs; i++)
    {
      freqVector[i] = Min_Freq + (Freq_Step * i);
    }


  srand((unsigned)time(0));
  
  cout << "Roulette Initialized" << endl;

  if(argc != 4)
    {cout << "Error: Usage. Specify start or cont, as well as NPop (Ex: Start 10)." << endl;}
  else
    {
      if(string(argv[1]) == "start")
	{
	  std::default_random_engine generator;
	  generator.seed(time(0));
	  // More gene things 
	  //
	  //
	  //

	  for(int i=0; i<NPop; i++)
	    {
	      for(int j=0; j< NSections; j++)
		{
		  for(int k=0; k<NVars; k++)
		    {
		      if(k == 0)
			{
			  // more genes references
			}
		    }
		}
	    }
	  dataWrite(Npop, varOutput, freq_coeffs, freqVector);
	  double meanTotal =0.0;
	  for(int i=0; i<Npop; i++)
	    {
	      meanTotal = meanTotal + varOutput[i][0][1];
	    }
	  float meanForGridSize = meanTotal / NPop;
	  ofstream datasize;
	  datasize.open("datasize.txt");
	  datasize << meanForGridSize/50.0 << ";";
	  datasize.close();
	}
      else if(string(argv[1]) == "cont")
	{
	  dataRead(varInput, fitness);
	  if (checkConvergence(varInput, fitness) == 1)
	    {
	      remove("highfive.txt");
	      ofstream highfive;
	      highfive.open("highfive.txt");
	      highfive << 1;
	      highfive.close();
	    }
	  else
	    {
	      roulette(varInput, varOutput, fitness);
	      cout << "Roulette complete." << endl;
	      dataWrite(NPop, varOutput, freq_coeff, freqVector);
	      double meanTotal =0.0;
	      for(int i=0; i<NPop; i++)
		{
		  meanTotal = meanTotal + varOutput[i][0][1];
		}
	      float meanForGridSize = meanTotal / NPop;
	      ofstream datasize;
	      datasize.open("datasize.txt");
	      datasize << meanForGridSize/50.0 << ";";
	      datasize.close();
	    }
	}
    }
  return (0);
}
