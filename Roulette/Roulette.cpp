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

void dataRead(vector<vector<vector<float>>> &varInput, vector<float> &fitness);

void dataWrite(int numChildren, vector<vector<vector<float>>> &varVector, int freq_coeffs, vector<double> freqVector);

void checkConvergence(vector<vector<vector<float>>> &varInput, vector<float> &fitness);

void roulette(vector<vector<vector<float>>> &varInput, vector<vector<vector<float>>> &varOutput, vector<float> &fitness);


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
//
//Gene 1 (controls angle 1 radians X axis)
float Initial_Mean_C1_G1 = M_PI ;
float Initial_Std_Dev_C1_G1 = M_PI / 2;

// Gene 2 (Y axis angle)
float Initial_Mean_C1_G2 = M_PI;
float Initial_Std_Dev_C1_G2 = M_PI / 2; 

// Gene 3 (Z axis)
float Initial_Mean_C1_G3 = M_PI;
float Initial_Std_Dev_C1_G3 =M_PI / 2;

// unused
const float Initial_Mean_CX_GY = 0.0f;
const float Initial_Std_Dev_CX_GY = 0.0f;
const float Mut_Modulator = 4.0f;


// Main Function

int main(int argc, char const *argv[])
{

  double Geoscale_Factor = stod(argv[3]);

  // frequencies scaled inversely with dimensions

  Max_Freq *= Geoscale_Factor;
  
  Min_Freq *= Geoscale_Factor;

  Freq_Step *= Geoscale_Factor;


  // Define NPop

  NPop = atoi(argv[2]);
 
  vector<vector<vector<float>>> varInput (NPop,vector<vector<float> >(NSections, vector <float>(NVars, 0.0f)));

  vector<float> fitness (NPop, 0.0f);

  vector<vector<vector<float>>> varOutput (NPop, vector<vector<float> >(NSections, vector <float>(NVars, 0.0f)));

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
	  generator.seed(1); // 1 is interchangable with time(0)
	  std::normal_distribution <float> distribution_X_Angle(Initial_Mean_C1_G1, Initial_Std_Dev_C1_G1);
	  std::normal_distribution <float> distribution_Y_Angle(Initial_Mean_C1_G2, Initial_Std_Dev_C1_G2);
	  std::normal_distribution <float> distribution_Z_Angle(Initial_Mean_C1_G3, Initial_Std_Dev_C1_G3);
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
			  float r = distribution_X_Angle(generator);
			  
			  while(r<=0)
			    r =distribution_X_Angle(generator);
			  varOutput[i][j][k] =r;
			}
		      else if (k==1)
			{
			  float r = distribution_Y_Angle(generator);
			  
			  while(r<=0)
			    r = distribution_Y_Angle(generator);
			  varOutput[i][j][k] = r;
			}
		      else if (k == 2)
			{
			  float r = distribution_Z_Angle(generator);

			  while(r <= 0)
			    r= distribution_Z_Angle(generator);
			  varOutput[i][j][k] = r;
			}
		    }
		}
	    }

	  /* std::default_random_engine generator;
	  generator.seed(time(0));
	  std::normal_distribution <float> distribution(Initial_Mean, Initial_Std_Dvn);
	  
	  for(int i=0; i<NPop; i++)
	    {
	      for(int j=0; i<NSections; j++)
		{
		  for(int k=0; k<NVars; k++)
		    {
		      float r = distribution(generator);
		      
		      while(r<=0)
			r = distribution(generator);
		      varOutput[i][j][k] = r;
		    }
		}
	    }
	  */

	  dataWrite(NPop, varOutput, freq_coeffs, freqVector);
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
	      dataWrite(NPop, varOutput, freq_coeffs, freqVector);
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
//
//
//
//
//
// Functions:

void dataWrite(int numChildren, vector<vector<vector<float>>> &varVector, int freq_coeffs, vector<double> freqVector)
{
  ofstream generationDNA;
  generationDNA.open("generationDNA.csv");
  generationDNA << "this is just a test";

  for(int i=0; i<freq_coeffs; i++)
    {
      if(i==freq_coeffs-1)
	{
	  generationDNA << freqVector[i] << "\n";
	}
      else
	{
	  generationDNA << freqVector[i] << ",";
	}
    }
  generationDNA << "Matricies for this generation: " << "\n";
  for(int i=0; i<numChildren; i++)
    {
      for(int j=0;j<NSections; j++)
	{
	  for(int k=0; k<NVars; k++)
	    {
	      if(k==(NVars-1))
		{
		  generationDNA << varVector[i][j][k] << "\n";
		}
		else
		  {
		    generationDNA << varVector[i][j][k] << ",";
		  }
	    }
	}
    }
      generationDNA.close();
}



void dataRead(vector<vector<vector<float>>> &varInput, vector<float> &fitness)
{
  ifstream generationDNA;
  generationDNA.open("generationDNA.csv");
  int csv_file_size = DNA_Garbage_End + (NPop * Nsections);
  string csvContent[csv_file_size+1];
  string strToDbl;

  for(int i=1; i<=csv_file_size; i++);
  {
    getline(generationDNA,csvContent[i]);
    if (i>DNA_Garbage_End)
      {
	double j=floor((i-DNA_Garbage_End)/NSections);
	int p=i-DNA_Garbage_End-NSections*j;
	istringstream stream(csvContent[i]);
	for(int k=0;k<NVars;k++)
	  {
	    getline(stream, strToDbl, ',');
	    varInput[j-1][p][k] = atof(strtoDbl.c_str());
	  }
      }
  }
  generationDNA.close();

  ifstream fitnessScores;
  fitnessScores.open("fitnessScores.csv");
  string fitnessRead[NPop+2];
  for(int i=0; i<(NPop+2); i++)
    {
      getline(fitnessScores, fitnessRead[i]);
      if(i>=2)
	{
	  fitness[i-2] = atof(fitnessRead[i].c_str());
	  if (fitness[i-2] < 0)
	    {
	      fitness[i-2] = 0;
	    }
	}
    }
  fitnessScores.close();
}


int checkConvergence(vector<vector<vector<float>>> &varInput, vector<float> &fitness)
{
  vector<vector<float>> meanTensor (NSections, vector<float>(NVars, 0));
  vector<vector<float>> dvnTensor (NSections, vector<float>(NVars,0));

  for(int i=0; i<Nvars; i++)
    {
      for(int j=0; j<NSections; j++)
	{
	  float totalSum=0.0f;
	  for (int k=0; k<NPop; k++)
	    {
	      totalSum+=varInput[k][j][i];
	    }
	  float mean = totalSum / NPop;
	  meanTensor[j][i] = mean;
	}
    }
  
  for(int i=0; i<NVars; i++)
    {
      for(int j=0; j<NSections; j++)
	{
	  float dvnDum = 0.0f;
	  for (int k=0; k<NPop;k++)
	    {
	      dvnSum+=pow((varInput[k][j][i] - meanTensor[j][i]),2);
	    }
	  float dvn = pow((dvnSum / (NPop -1)), 1/2);
	  dvnTensor[j][i] = dvn;
	}
    }
  if (dvnTensor[0][1] <= Convergence)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}


void roulette(vector<vector<vector<float>>> &varInput, vector<vector<vector<float>>> &varOutput, vector<float> &fitness)
{
  float fitness_total =0;
  int roulette_no=NPop*(1-Tourney_Proportion);
  int tourney_no= NPop*Tourney_Proportion;

  while (roulette_no+tourney_no!=NPop)
    {
      if(roulette_no+tourney_no>NPop && roulette_no>tourney_no)
	{
	  tourney_no+=-1;
	}
      else if (roulette_no+tourney_no <NPop && roulette_no>tourney_no)
	{
	  roulette_no+=1;
	}
      else if (roulette_no+tourney_no<NPop && roulette_no<tourney_no)
	{
	  tourney_no+=1;
	}
    }

  vector<vector<int>> selected (roulette_no, vector<int>(Parent_no ,0));
  
  for(int i=0; i<NPop; i++)
    {
      fitness_total+=fitness[i];
    }
  
  for(int i=0; i<roulette_no;i++)
    {
      for(int j=0; j<Parent_no; j++)
	{
	  float partial_sum =0.0f;
	  floar r =rand()%100;
	  
	  for(int k=0; k<NPop; k++)
	    {
	      partial_sum = partial_sum + fitness[k];
	      if(r<(100*partial_sum/fitness_total))
		{
		  selected[i][j]=k;
		  break;
		}
	    }
	}
    }




  for(int i=0; i<roulette_no; i++)
    {
      for(int j=0; j<NSections; j++)
	{
	  for(int k=0; k<NVars;k++)
	    {
	      int pick1=rand()%Parent_no;
	      varOutput[i][j][k] = varInput[selected[i][pick1]][j][k];
	    }
	}
    }

  vector<vector<float>> meanTensor (NSections, vector<float>(NVars,0));
  vector<vector<float>> dvnTensor (NSections, vector<float>(NVars,0));

  for (int i=0; i<NVars; i++)
    {
      for( int j=0; j<NSections; j++)
	{
	  float totalSum=0.0f;
	  for (int k=0; k<NPop; k++)
	    {
	      totalSum +=varInput[k][j][i];
	    }
	  float mean = totalSum / NPop;
	  meanTensor[j][i]=mean;
	}
    }

  for (int i=0; i<NVars; i++)
    {
      for (int j=0; j<NSections; j++)
	{
	  float dvnSum = 0.0f;
	  for (int k=0; k<NPop;k++)
	    {
	      dvnSum+=pow((varInput[k][j][i] - meanTensor[j][i]),2);
	    }
	  float dvn = pow((dvnSum / (NPop - 1)),1/2);
	  dvnTensor[j][i]=dvn;
	}
    }


  vector<bool> mutate_flag (roulette_no, false);

  int roul_mut = roulette_no *Mutability;

  for(int i=0; i<roul_mut; i++)
    {
      int r =rand()%roulette_no;

      while(mutate_flag[r] == true)
	{
	  r = rand()%roulete_no;
	}
      for (int j=0; j<NSections; j++)
	{
	  int numberOfMutations = rand()%Nvars +1;

	  for(int k=0; k<numberOfMutations; k++)
	    {
	      int chromosomeMutation = rand%NSections;
	      int geneMutation = rand()%NVars;
	      std::default_random_engine generator;
	      generator.seed(time(0));
	      std::normal_distribution <float> distribution(meanTensor[chromosomeMutation][geneMutation], dvnTensor[chromosomeMutation][geneMutation]);
	      
	      int coefficient=rand()%2;
	      if(coefficient == 0)
		{
		  varOutput[r][chromosomeMutation][geneMutation]=varOutput[r][chromosomeMutation][geneMutation]+(distribution(generator)/Mut_Modulator);
		}
	      else
		{
		  varOutput[r][chromosomeMutation][geneMutation]=varOutput[r][chromosomeMutation][geneMutation]-(distribution(generator)/Mut_Modulator);
		}
	      while(varOutput[r][chromosomeMutation][geneMutation]<=0)
		{
		  varOutput[r][chromosomeMutation][geneMutation]=varOutput[r][chromosomeMutation][geneMutation]+distribution(generator);
		}
		mutate_flag[r] =true;
	    }
	}
    }
}

