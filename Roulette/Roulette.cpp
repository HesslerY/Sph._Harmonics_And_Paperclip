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
