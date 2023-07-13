#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "mpi.h"
#include "Random/random.h"
#include "Population.hpp"
#include "Individual.hpp"
#include "Position2D.hpp"
#include "utils.hpp"

using namespace std;



int main(int argc, char** argv){

    // initialize mpi
    int size, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    

    Random random(rank + 1);

    int nGens, migrationFreq, nIndividuals;
    bool exchange;

    ifstream input("input.in");
    input >> exchange >> nGens >> migrationFreq >> nIndividuals;
    input.close();

    vector<Position2D> positions = readCapitals();

    int nCities = positions.size();

    // in this case pre-computing distances actually speeds up the code
    // by about 25% (this was found with some VERY rough estimates)
    // since we use trigonometric functions to compute distances
    vector<double> distances = computeDistances(nCities, 1., positions);

    Population population(nIndividuals, nCities, random, positions, distances);

    string path = "independent/";
    if (exchange) path = "exchange/";


    for (int i = 0; i < nGens; i++){
        population.nextGeneration(random);

        // the best score/path is evaluated in the previous function before crossover and mutations
        // these following print statements refer to the population before calling nextGeneration()
        population.printBestPath(path + "path_" + to_string(rank) + ".out");
        population.printAveragePath(path + "average_" + to_string(rank) + ".out");

        // we want to pass the best element of each process to another random one
        // to do this we exploit what we already have: we generate an "Individual" object
        // of length=("size"+1) randomly ("size" is the number of processes)
        // since the first element is fixed to 1 in an "Individual", we remove the first
        // element and subtract 2 to all the other ones, leaving us with a vector of integers 
        // from 0 to "size"-1 called "sendTo"
        // then we pass the best element of process i to process "sendTo[i]"
        // and process j will receive the best element from "receiveFrom[j]"
        if (exchange){
            if ((i + 1) % migrationFreq == 0 && i + 1 < nGens){
                // decides which processes exchange information as stated above (only once)
                vector<int> sendTo(size, -1);
                vector<int> receiveFrom(size, -1);

                if (rank == 0){
                    Individual idk(size + 1, random);
                    sendTo = idk.getGenes();

                    sendTo.erase(sendTo.begin());
                    for (int i = 0; i < size; i++) sendTo[i] -= 2;

                    for (int i = 0; i < size; i++) receiveFrom[sendTo[i]] = i;
                }

                MPI_Barrier(MPI_COMM_WORLD); // wait until "sendTo" and "receiveFrom" are generated
                MPI_Bcast(&sendTo[0], size, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&receiveFrom[0], size, MPI_INT, 0, MPI_COMM_WORLD);
                

                population.computeScores();
                population.order();

                vector<int> bestToMigrate = population.getGenes(0);
                int* newGenes = new int[nCities];

                int destination = sendTo[rank];
                int source = receiveFrom[rank];
                
                MPI_Status status;
                MPI_Request request;
                
                MPI_Isend(&bestToMigrate[0], nCities, MPI_INT, destination, rank, MPI_COMM_WORLD, &request);
                MPI_Recv(newGenes, nCities, MPI_INT, source, source, MPI_COMM_WORLD, &status);
                
                // before overwriting the best element, wait until every process receives "newGenes"
                // otherwise it might send the genes it received from another process and not its own best.
                MPI_Barrier(MPI_COMM_WORLD);
                population.setIndividual(0, Individual(vector<int>(newGenes, newGenes + nCities)));

                if (rank == 0) cout << "generation " << i + 1 << ": migration completed" << endl;
            }
        }
    }

    // sort and print the last generation
    population.computeScores();
    population.order();
    population.printBestPath(path + "path_" + to_string(rank) + ".out");
    population.printAveragePath(path + "average_" + to_string(rank) + ".out");


    cout << "process " << rank << " done\n";
    population.print(0, path + "best.out");



    MPI_Finalize();

    return 0;
}