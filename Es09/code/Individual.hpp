#ifndef __Individual__
#define __Individual__

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <iomanip>
#include "Random/random.h"
#include "utils.hpp"
#include "Position2D.hpp"

using namespace std;


class Individual{

public:
    Individual(){genes = {}; size = 0; score = 0.;}
    Individual(vector<int> g){genes = g; size = genes.size(); score = 0.;}
    Individual(int, Random&);
    Individual(const Individual& i){genes = i.genes; size = genes.size(); check(); score = i.score;}
    ~Individual(){}

    Individual operator=(const Individual&);

    int getSize() const {return size;}
    vector<int> getSlice(int start, int end) const {return vector<int>(genes.begin() + start, genes.begin() + end);}
    vector<int> getGenes() const {return genes;}
    // score is determined externally since the methods can vary
    // but it is useful to have it assigned directly to each individual
    void setScore(double s){score = s;}
    double getScore() const {return score;}

    // operators
    void swapPair(int, int);
    void swapPair(Random&);
    void swapNeighbours(Random&);
    void flipBlock(int, int);
    void flipBlock(Random&);
    void swapBlock(int, int, int);
    void swapBlock(Random&);
    void moveBlock(int, int, int);
    void moveBlock(Random&);
    
    int find(int);
    void check();
    void print();
    void print(string);

private:
    vector<int> genes;
    int size;
    double score;
};

#endif