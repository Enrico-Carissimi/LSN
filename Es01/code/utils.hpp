#ifndef __utils__
#define __utils__

#include "Random/random.h"

#include <cmath>
#include <vector>

using namespace std;

double calculateError(double sum, double sum2, int N){
    return (N > 1) ? sqrt(((sum2 - sum * sum / N) / N) / (N - 1)) : 0;
}

double chi2(vector<int> O, double E){ // (in this case E is always the same, E = 100)
    double sum = 0;

    for (size_t i = 0; i < O.size(); i++){
        sum += (O[i] - E) * (O[i] - E) / E;
    }

    return sum;
}

double meanUniform(int N, Random& random){
    double sum = 0;

    for (int i = 0; i < N; i++){
        sum += random.Rannyu();
    }

    return sum / N;
}

double meanExp(int N, Random& random){
    double sum = 0;

    for (int i = 0; i < N; i++){
        sum += random.Exp(1.);
    }

    return sum / N;
}

double meanLorentz(int N, Random& random){
    double sum = 0;

    for (int i = 0; i < N; i++){
        sum += random.Lorentz(0., 1.);
    }

    return sum / N;
}

double randomSineWOPi(Random& random){
    double x, y, r2;
    do{
        x = random.Rannyu(-1., 1.);
        y = random.Rannyu(-1., 1.);
        r2 = x * x + y * y;
    } while (r2 > 1.);

    return y / sqrt(r2);
}

#endif