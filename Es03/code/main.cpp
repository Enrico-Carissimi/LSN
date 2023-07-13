#include "Random/random.h"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int main(){
    double const S0 = 100., T = 1., K = 100., r = 0.1, sigma = 0.25; // S0 is the initial price of the asset

    Random random;

    int numberThrows = 100000;
    int numberBlocks = 100;
    int throwsPerBlock = numberThrows / numberBlocks; // numberThrows should be a multiple of numberBlocks

//es 1.1

    ofstream dataOutput1("data/data-3.1.1-call.out");
    ofstream dataOutput2("data/data-3.1.1-put.out");

    double sumOverBlocksCall = 0, sumSquaredOverBlocksCall = 0;
    double sumOverBlocksPut = 0, sumSquaredOverBlocksPut = 0;

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double sumCall = 0., sumPut = 0.;

        for (int i = 0; i < throwsPerBlock; i++){
            // we use a discrete GBM with one step: t_0 = 0 and t_1 = T
            // to sample directly the final price S(T)
            // we can multiply by the discount e^(-rT) at the end since here r and T are constants
            sumCall += call(discreteGBM(1, S0, r, sigma, T, random), K);
            sumPut += put(discreteGBM(1, S0, r, sigma, T, random), K);
        }
        sumCall /= throwsPerBlock, sumPut /= throwsPerBlock;

        sumOverBlocksCall += sumCall;
        sumSquaredOverBlocksCall += sumCall * sumCall;
        sumOverBlocksPut += sumPut;
        sumSquaredOverBlocksPut += sumPut * sumPut;

        dataOutput1 << sumOverBlocksCall * exp(-r * T) / blockIndex << "\t" << calculateError(sumOverBlocksCall, sumSquaredOverBlocksCall, blockIndex) << "\n";
        dataOutput2 << sumOverBlocksPut * exp(-r * T) / blockIndex << "\t" << calculateError(sumOverBlocksPut, sumSquaredOverBlocksPut, blockIndex) << "\n";
    }

    dataOutput1.close();
    dataOutput2.close();

//es 1.2

    dataOutput1.open("data/data-3.1.2-call.out");
    dataOutput2.open("data/data-3.1.2-put.out");

    sumOverBlocksCall = 0, sumSquaredOverBlocksCall = 0;
    sumOverBlocksPut = 0, sumSquaredOverBlocksPut = 0;

    int numberIntervals = 100;
    double dt = T / numberIntervals;

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double sumCall = 0., sumPut = 0.;

        for (int i = 0; i < throwsPerBlock; i++){
            // we can multiply by the discount e^(-rT) at the end since here r and T are constants
            sumCall += call(discreteGBM(numberIntervals, S0, r, sigma, dt, random), K);
            sumPut += put(discreteGBM(numberIntervals, S0, r, sigma, dt, random), K);
        }
        sumCall /= throwsPerBlock, sumPut /= throwsPerBlock;

        sumOverBlocksCall += sumCall;
        sumSquaredOverBlocksCall += sumCall * sumCall;
        sumOverBlocksPut += sumPut;
        sumSquaredOverBlocksPut += sumPut * sumPut;

        dataOutput1 << sumOverBlocksCall * exp(-r * T) / blockIndex << "\t" << calculateError(sumOverBlocksCall, sumSquaredOverBlocksCall, blockIndex) << "\n";
        dataOutput2 << sumOverBlocksPut * exp(-r * T) / blockIndex << "\t" << calculateError(sumOverBlocksPut, sumSquaredOverBlocksPut, blockIndex) << "\n";
    }

    dataOutput1.close();
    dataOutput2.close();

    return 0;
}