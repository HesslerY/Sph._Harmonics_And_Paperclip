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

int main(){
	double array[] = {.1, 0, 45, .9, .7536};
	
	for (int i = 0; i <= 4; i++){
		cout << array[i] << endl;
	}
	
	insertionSort(array, 5);
	
	for (int i = 0; i <= 4; i++){
		cout << array[i] << endl;
	}
}
