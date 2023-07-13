#include <iostream>
#include <fstream>
#include "utils.hpp"

using namespace std;

double potential(double x){
    return pow(x, 4) - 2.5 * x * x;
}



double psiT(double x, double sigma, double mu){
    double d = 0.5 / (sigma * sigma);
    double n = (x - mu);
    double N = (x + mu);

    return exp(- n * n * d) + exp(- N * N * d);
}



double d2psiT(double x, double sigma, double mu){
    double a = pow((x - mu)/sigma, 2);
    double A = pow((x + mu)/sigma, 2);

    return (1. / (sigma * sigma)) * (-psiT(x, sigma, mu) + a * exp(- a * 0.5) + A * exp(- A * 0.5));
}



double error(double sum, double sum2, int N){
    return (N > 1) ? sqrt(((sum2 - sum * sum / N) / N) / (N - 1)) : 0;
}



// generates a random point within sigma of one of the two centres (mu or -mu)
double generateStartingPoint(double sigma, double mu, Random& random){
    double x, p;
    p = (random.Rannyu() > 0.5) ? -mu : mu;
    x = (p + (random.Rannyu() - 0.5) * 2. * sigma);
    return x;
}