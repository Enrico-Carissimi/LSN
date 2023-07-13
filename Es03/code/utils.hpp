#ifndef __utils__
#define __utils__

#include "Random/random.h"

#include <cmath>
#include <algorithm>

using namespace std;

double calculateError(double sum, double sum2, int N){
    return (N > 1) ? sqrt(((sum2 - sum * sum / N) / N) / (N - 1)) : 0;
}

double discreteGBM(int nstep, double S0, double mu, double sigma, double dt, Random& random){
    double S = S0;

    for (int i = 0; i < nstep; i++){
        S = S * exp((mu - 0.5 * sigma * sigma) * dt + sigma * random.Gauss(0., 1.) * sqrt(dt));
    }

    return S;
}

double call(double S, double K){
    return max(0., S - K);
}

double put(double S, double K){
    return max(0., K - S);
}

#endif