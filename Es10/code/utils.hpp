#ifndef __utils__
#define __utils__

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include "Random/random.h"
#include "Position2D.hpp"

using namespace std;

vector<Position2D> readCapitals();

template <typename T> double distance2(const T& a, const T& b){return (a - b) * (a - b);}
double distanceOnSphere(double r, Position2D, Position2D);
vector<double> computeDistances(int, double r, vector<Position2D>);

void printPositions(string, vector<Position2D>);
template <typename T> void printVector(vector<T> v){for (int i = 0; i < (int)v.size(); i++) cout << v[i] << " "; cout << endl;}

#endif