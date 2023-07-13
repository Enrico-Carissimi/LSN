#include "Random/random.h"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int main(){
    Random random;

//es 1.1

    int numberThrows = 100000;
    int numberBlocks = 100;
    int throwsPerBlock = numberThrows / numberBlocks; //numberThrows should be a multiple of numberBlocks
    
    ofstream dataOutput("data/data-1.1.1.out");

    double sumOverBlocks = 0, sumSquaredOverBlocks = 0;

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double sum = 0;

        for (int i = 0; i < throwsPerBlock; i++){
            sum += random.Rannyu();
        }
        sum /= throwsPerBlock;

        sumOverBlocks += sum;
        sumSquaredOverBlocks += sum * sum;

        dataOutput << sumOverBlocks / blockIndex << "\t" << calculateError(sumOverBlocks, sumSquaredOverBlocks, blockIndex) << "\n";
    }

    dataOutput.close();

//es 1.2

    dataOutput.open("data/data-1.1.2.out");

    sumOverBlocks = 0, sumSquaredOverBlocks = 0;

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double sum = 0;

        for (int i = 0; i < throwsPerBlock; i++){
            sum += pow(random.Rannyu() - 0.5, 2);
        }
        sum /= throwsPerBlock;

        sumOverBlocks += sum;
        sumSquaredOverBlocks += sum * sum;

        dataOutput << sumOverBlocks / blockIndex << "\t" << calculateError(sumOverBlocks, sumSquaredOverBlocks, blockIndex) << "\n";
    }

    dataOutput.close();

//es 1.3

    int M = 100; // number of intervals
    int n = 10000; // number of throws per interval
    int numberRuns = 100;
    double expectedOccurences = (double) n / M;

    dataOutput.open("data/data-1.1.3.out");

    for (int j = 0; j < numberRuns; j++){
        vector<int> counts(M, 0);

        for (int i = 0; i < n; i++){
            double r = random.Rannyu();
            counts[int(r * M)]++; // r*M is truncated, so that [0, 0.01) -> 0, [0.01, 0.02) -> 1, ... , [0.99, 1) -> 99
        }

        dataOutput << chi2(counts, expectedOccurences) << "\n";
    }

    dataOutput.close();

//es 2

    numberRuns = 10000;

    dataOutput.open("data/data-1.2-uniform.out");

    for (int i = 0; i < numberRuns; i++){
        dataOutput << meanUniform(1, random) << "\t" << meanUniform(2, random) << "\t" << meanUniform(10, random) << "\t" << meanUniform(100, random) << "\n";
    }

    dataOutput.close();

    dataOutput.open("data/data-1.2-exp.out");

    for (int i = 0; i < numberRuns; i++){
        dataOutput << meanExp(1, random) << "\t" << meanExp(2, random) << "\t" << meanExp(10, random) << "\t" << meanExp(100, random) << "\n";
    }

    dataOutput.close();

    dataOutput.open("data/data-1.2-lorentz.out");

    for (int i = 0; i < numberRuns; i++){
        dataOutput << meanLorentz(1, random) << "\t" << meanLorentz(2, random) << "\t" << meanLorentz(10, random) << "\t" << meanLorentz(100, random) << "\n";
    }

    dataOutput.close();

//es 3

    //we consider 2 parallel lines in 2d at y = 0 and y = d (the second is actually excluded, to have an infinite pattern)
    //the needle is simulated as 2 point at distance L
    //we might generate 2 pseudo random numbers: the y position of the first point (y1, in the interval [0, d))
    //and the angle between the needle and the x axis (in the interval [0, 2pi))
    //but if we extract the angle directly, this method actually uses pi to calculate pi

    //we are interested in the y position of the second end of the needle, y1 + L*sin(theta), where theta is the angle discussed above
    //so we generate directly the sine:
    //we can generate a point inside a circle of unit radius and then normalize its distance from the centre to 1

    //if this segment (the needle) intersects the lines, then we have a positive outcome
    //so we count how many times the second point is outside the [0, d) interval

    double d = 10.;
    double L = 5.;

    dataOutput.open("data/data-1.3.out");

    numberThrows = 10000000, numberBlocks = 250;
    throwsPerBlock = numberThrows / numberBlocks;

    sumOverBlocks = 0, sumSquaredOverBlocks = 0;

    for (int blockIndex = 1; blockIndex <= numberBlocks; blockIndex++){
        double numberIntersections = 0;

        for (int i = 0; i < throwsPerBlock; i++){
            double y1 = random.Rannyu(0., d);

            double y2 = y1 + L * randomSineWOPi(random);

            if (y2 >= d || y2 <= 0) numberIntersections++;
        }
        
        double pi = 2. * L * throwsPerBlock / (d * numberIntersections);

        sumOverBlocks += pi;
        sumSquaredOverBlocks += pi * pi;


        dataOutput << sumOverBlocks / blockIndex << "\t" << calculateError(sumOverBlocks, sumSquaredOverBlocks, blockIndex) << "\n";
    }

    dataOutput.close();

    return 0;
}