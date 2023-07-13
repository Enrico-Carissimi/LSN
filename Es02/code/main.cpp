#include "Random/random.h"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int main(){
    int numberThrows = 10000, numberBlocks = 100;
    int throwsPerBlock = numberThrows / numberBlocks; // numberThrows must be a multiple of numberBlocks!

    Random random;

    ofstream dataOutput("data/data-2.1.1.out");
        
    double sumOverBlocks = 0., sumSquaredOverBlocks = 0.;

//es 1.1

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double sum = 0.;

        for (int i = 0; i < throwsPerBlock; i++){
            sum += integrand(random.Rannyu());
        }

        sum /= throwsPerBlock;

        sumOverBlocks += sum, sumSquaredOverBlocks += sum * sum;

        dataOutput << sumOverBlocks / blockIndex << "\t" << calculateError(sumOverBlocks, sumSquaredOverBlocks, blockIndex) << "\n";
    }

    dataOutput.close();

//es 1.2

    // we use 2-2x as the probability distribution and we integrate integrand/(2-2x)
    // the cumulative function is: y = 2x-x^2
    // thus we should sample x = 1+sqrt(1-y)

    dataOutput.open("data/data-2.1.2.out");

    sumOverBlocks = 0, sumSquaredOverBlocks = 0;

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double sum = 0;

        for (int i = 0; i < throwsPerBlock; i++){
            double x = probabilityDistribution(random);
            sum += integrand(x) * 0.5 / (1 - x);
        }

        sum /= throwsPerBlock;

        sumOverBlocks += sum, sumSquaredOverBlocks += sum * sum;

        dataOutput << sumOverBlocks / blockIndex << "\t" << calculateError(sumOverBlocks, sumSquaredOverBlocks, blockIndex) << "\n";
    }

    dataOutput.close();

//es 2.1 (la costante k è 1)

    int numberSteps = 100;

    int numberRuns = 10000;
    numberBlocks = 100;
    int runsPerBlock = numberRuns / numberBlocks; // numberRuns must be a multiple of numberBlocks!

    double a = 1.;

    vector<double> sumOverSteps(numberSteps + 1, 0.);
    vector<double> sumSquaredOverSteps(numberSteps + 1, 0.);
    Position origin; //initialized automatically to (0, 0, 0)

    dataOutput.open("data/data-2.2.1.out");

    // we simulate numberRuns random walks on a cubic lattice
    // for each one, at every step we calculate the distance from the origin and accumulate it
    // we then average over the distances at each step (using block averages)
    // we print only the final value with its error
    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        vector<double> sum(numberSteps + 1, 0.);

        for (int j = 0; j < runsPerBlock; j++){
            Position pos;

            for (int step = 0; step <= numberSteps; step++){
                sum[step] += pos.distance2(origin);
                pos.move(moveCubicLattice(a, random));
            }
        }

        for (int step = 0; step <= numberSteps; step++){
            double d = sum[step] / runsPerBlock;

            sumOverSteps[step] += sqrt(d);
            sumSquaredOverSteps[step] += d; // = sqrt(d)^2
        }
    }

    for (int step = 0; step <= numberSteps; step++){
        dataOutput << sumOverSteps[step] / numberBlocks << "\t" << calculateError(sumOverSteps[step], sumSquaredOverSteps[step], numberBlocks) << "\n";
    }

    dataOutput.close();

//es 2.2

    // reset
    fill(sumOverSteps.begin(), sumOverSteps.end(), 0.);
    fill(sumSquaredOverSteps.begin(), sumSquaredOverSteps.end(), 0.);

    dataOutput.open("data/data-2.2.2.out");

    // we do the same as before, but in the continuum
    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        vector<double> sum(numberSteps + 1, 0.);

        for (int j = 0; j < runsPerBlock; j++){
            Position pos;

            for (int step = 0; step <= numberSteps; step++){
                sum[step] += pos.distance2(origin);
                pos.move(moveContinuum(a, random));
            }
        }

        for (int step = 0; step <= numberSteps; step++){
            double d = sum[step] / runsPerBlock;

            sumOverSteps[step] += sqrt(d);
            sumSquaredOverSteps[step] += d; // = sqrt(d)^2
        }
    }

    for (int step = 0; step <= numberSteps; step++){
        dataOutput << sumOverSteps[step] / numberBlocks << "\t" << calculateError(sumOverSteps[step], sumSquaredOverSteps[step], numberBlocks) << "\n";
    }

    dataOutput.close();

    return 0;
}