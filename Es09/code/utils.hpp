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

vector<double> generateRandom(int, double min, double max, Random&);
vector<Position2D> generateRandom2D(int, double r, Random&); // on circle of radius r
vector<Position2D> generateRandom2D(int, double xmin, double xmax, double ymin, double ymax, Random&); //in a rectangle

vector<Position2D> generateOnCircle(vector<double>);

template <typename T> double distance2(const T& a, const T& b){return (a - b) * (a - b);}
double distanceOnCircle(double r, double, double);

vector<double> computeDistancesCircle(int, double, vector<double>);
vector<double> computeDistancesSquare(int, vector<Position2D>);

void printPositions(string, vector<Position2D>);
template <typename T> void printVector(vector<T> v){for (int i = 0; i < (int)v.size(); i++) cout << v[i] << " "; cout << endl;}

#endif