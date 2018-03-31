#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <vector>
using namespace std;

double Pi = 3.14159265359;
int PopMAX = 100;


int main()
{
	
int Sph[13];
int SphHarMax = 6;

for (int i = 0; i < 13; i++){
	if (i < SphHarMax){
		Sph[i]=1;
	}
	else{
		Sph[i]=0;
	}
}

for (int i = 0; i < 13; i++){
	cout << i << " : " << Sph[i] << endl;
}

return 0;
}
