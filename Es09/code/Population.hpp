#ifndef __Population__
#define __Population__

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include "utils.hpp"
#include "Individual.hpp"
#include "Position2D.hpp"
#include "Random/random.h"

using namespace std;



class Population{

public:
    Population(vector<Individual>, vector<Position2D>, vector<double>);
    Population(int, int, Random&, vector<Position2D>, vector<double>);
    ~Population(){}

    int getSize(){return size;}
    int getIndividualSize(){return individualSize;}
    double getScore(int i){return population[i].getScore();}
    double getBestScore(){return bestScore;}

    void add(Individual i){population.push_back(i); size++; check();}
    void check();
    void setIndividualScore(int);
    void computeScores();
    void computeAveragePath();
    void order();

    void mutateIndividual(int, Random&);
    void mutate(Random&);
    Individual* crossover(int, int, int);
    int selection(Random&);
    void nextGeneration(Random&);

    void print(int i){population[i].print();}
    void print(int i, string file){population[i].print(file);}
    void print(string file){for (int i = 0; i < size; i++) population[i].print(file);}
    void printBest(){bestIndividual.print();}
    void printBestPath(string file){ofstream out(file, ios::app); out << 1. / bestScore << endl; out.close();}
    void printAveragePath(string file){ofstream out(file, ios::app); out << averagePath << endl; out.close();}

private:
    vector<Individual> population;
    int size, individualSize;
    vector<Position2D> positions;
    vector<double> distances;
    double bestScore, averagePath;
    Individual bestIndividual;
};

#endif