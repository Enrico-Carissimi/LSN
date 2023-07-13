#include "Population.hpp"


Population::Population(vector<Individual> v, vector<Position2D> pos, vector<double> dist){
    population = v;
    size = population.size();
    individualSize = population[0].getSize();
    positions = pos;
    distances = dist;
    
    computeScores();

    check();
}



Population::Population(int popSize, int l, Random& random, vector<Position2D> pos, vector<double> dist){
    individualSize = l;
    size = popSize;
    positions = pos;
    distances = dist;

    vector<Individual> v;
    for (int i = 0; i < size; i++){
        v.push_back(Individual(individualSize, random));
    }
    population = v;

    computeScores();

    check();
}



void Population::check(){
    for (int i = 0; i < size; i++){

        if (population[i].getSize() != individualSize){
            cout << "Illegal individual: wrong size (" << population[i].getSize() << ", is " << individualSize << ")" << endl;
            cout << "\t";
            population[i].print();
            exit(-1);
        }

        population[i].check();
    }
}



void Population::setIndividualScore(int j){
    double sum = 0.;

    vector<int> v = population[j].getGenes();
    v.push_back(1);

    for (int i = 0; i < individualSize; i++){
        sum += distances[(v[i] - 1) * individualSize + v[i + 1] - 1];
    }

    population[j].setScore(1. / sum); // score is the inverse so it is higher for better solutions
}



void Population::computeScores(){
    for (int i = 0; i < size; i++){
        setIndividualScore(i);
    }
}



void Population::computeAveragePath(){
    double sum = 0;
    int n = size / 2;
    
    for (int i = 0; i < n; i++){
        sum += 1. / getScore(i);
    }

    averagePath = sum / n;
}



void Population::order(){
    // uses std functions to sort the array based on individual scores
    sort(population.begin(), population.end(), [](Individual a, Individual b){return a.getScore() < b.getScore();});
    reverse(population.begin(), population.end());

    bestScore = getScore(0);
    bestIndividual = population[0];
}



// declared in Population since it might depend on probabilities external to the individual
void Population::mutateIndividual(int i, Random& random){
        if (random.Rannyu() < 0.05) population[i].swapBlock(random);
        if (random.Rannyu() < 0.07) population[i].moveBlock(random);
        if (random.Rannyu() < 0.07) population[i].flipBlock(random);
        if (random.Rannyu() < 0.1) population[i].swapPair(random);
        if (random.Rannyu() < 0.12) population[i].swapNeighbours(random);
}



void Population::mutate(Random& random){
    for (int i = 0; i < size; i++){
        mutateIndividual(i, random);
    }
}



Individual* Population::crossover(int i, int j, int cutPosition){
    Individual* siblings = new Individual[2];

    if (cutPosition <= 0 || cutPosition >= individualSize){
        cout << "Illegal action: trying to cut individuals at out of bounds location" << endl;
        exit(-1);
    }

    vector<int> body1 = population[i].getSlice(0, cutPosition);
    vector<int> body2 = population[j].getSlice(0, cutPosition);
    vector<int> tail1 = population[i].getSlice(cutPosition, individualSize);
    vector<int> tail2 = population[j].getSlice(cutPosition, individualSize);

    vector<int> buffer1(individualSize, -1);
    vector<int> buffer2(individualSize, -1);

    for (int k = 0; k < individualSize - cutPosition; k++){
        int i1 = population[j].find(tail1[k]);
        int i2 = population[i].find(tail2[k]);

        buffer1[i1] = tail1[k];
        buffer2[i2] = tail2[k];
    }

    for (int k = individualSize - 1; k >= 0; k--){
        if (buffer1[k] <= 0) buffer1.erase(buffer1.begin() + k);
        if (buffer2[k] <= 0) buffer2.erase(buffer2.begin() + k);
    }

    body1.insert(body1.end(), buffer1.begin(), buffer1.end());
    body2.insert(body2.end(), buffer2.begin(), buffer2.end());

    siblings[0] = Individual(body1);
    siblings[1] = Individual(body2);

    return siblings;
}



int Population::selection(Random& random){
    return (int)(size * pow(random.Rannyu(), 4));
}



void Population::nextGeneration(Random& random){
    computeScores();
    order();

    bestIndividual = population[0];
    bestScore = getScore(0);
    computeAveragePath();

    vector<Individual> nextGen = {};

    for (int i = 0; i < size / 2; i++){
        if (random.Rannyu() < 0.75){
            int cutLocation = random.Integer(1, individualSize);
            Individual* siblings = crossover(selection(random), selection(random), cutLocation);
            nextGen.push_back(siblings[0]);
            nextGen.push_back(siblings[1]);
        }
    }

    int missing = size - nextGen.size();
    for (int i = 0; i < missing; i++){
        nextGen.push_back(population[selection(random)]);
    }

    population = nextGen;
    
    mutate(random);
}