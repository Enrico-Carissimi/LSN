#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include "Random/random.h"
#include "Population.hpp"
#include "Individual.hpp"
#include "Position2D.hpp"
#include "utils.hpp"

using namespace std;



int main(){
    Random random;

    int config, mode, nCities, nGens, nIndividuals;

    ifstream input("input.in");
    input >> config >> mode >> nCities >> nGens >> nIndividuals;
    input.close();

    string path;
    
    vector<Position2D> cities;

    // we pre-compute distances so we don't have to create two different "population" objects
    // that compute them in two different ways; this way we can use the same variables for both configurations
    // also, in this way we compute only N^2 distances (N=number of cities)
    // if we compute distances inside the population class, we would need,
    // for G generations and C cromosomes (size of the population), N*C*G > N^2 calculations
    // I don't think this makes actually any difference in the speed of this code, since we still
    // need to access the array of distances that many times, and in this case distances contain
    // only sums and products which are very fast to compute

    vector<double> distances;

    if (!config){ // circle
        path = "circle/";
        vector<double> angles = generateRandom(nCities, 0., 2 * M_PI, random);
        cities = generateOnCircle(angles);
        distances = computeDistancesCircle(nCities, 1., angles);
    }
    else{ //square
        path = "square/";
        cities = generateRandom2D(nCities, 0., 1., 0., 1., random);
        distances = computeDistancesSquare(nCities, cities);
    }

    printPositions(path + "positions.out", cities);

    Population population(nIndividuals, nCities, random, cities, distances);

    for (int i = 0; i < nGens; i++){
        population.nextGeneration(random);

        // the best score/path is evaluated in the previous function before crossover and mutations
        // these following print statements refer to the population before calling nextGeneration()
        population.printBestPath(path + "path.out");
        population.printAveragePath(path + "average.out");

        if (mode){
            cout << "generation " << i + 1 << " done\n";
            cout << "the best score before mutation was " << population.getBestScore();
            cout << "\ncorresponding to:\n";
            population.printBest();
            cout << endl;
        }
    }

    // sort and print the last generation
    population.computeScores();
    population.order();
    population.printBestPath(path + "path.out");
    population.printAveragePath(path + "average.out");


    if (mode) cout << "\n\n";

    cout << "the best candidate is:\n";
    population.print(0);
    population.print(0, path + "best.out");
    cout << "with score " << population.getBestScore() << "\nand path length of " << 1. / population.getBestScore() << endl;

    cout << "\nprinted best path for each generation to file \"" << path << "path.out\"" << endl;
    cout << "printed average path of the best half of the population for each generation to file \"" << path << "average.out\"" << endl;

    if (mode){
        population.print(path + "population.final");
        cout << "printed final population to file \"" << path << "population.final\"" << endl;
    }


    return 0;
}