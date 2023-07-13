#ifndef __utils__
#define __utils__

#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include "Random/random.h"

using namespace std;

double potential(double x);
double psiT(double x, double sigma, double mu);
double d2psiT(double x, double sigma, double mu);

double error(double sum, double sum2, int N);

double generateStartingPoint(double sigma, double mu, Random& random);



template <typename T>
vector<T> loadFile(string name){
    ifstream in(name);
    vector<T> data;
    T x;

    if (!in){
        cout << "Impossible to open file \"" << name << "\"" << endl;
        exit(-1);
    }

    while (!in.eof()){
        in >> x;
        if (!in.eof()) data.push_back(x);
    }

    in.close();

    return data;
}



template<typename Data, typename Func1, typename Func2>
bool metropolisStep(Data& x, Func1 p, Func2 T, Random& random){
    double pOld = p(x);
    double xOld = x;
    x = T(x);
    double pNew = p(x);

    if (random.Rannyu() <= min(1., pNew / pOld)){
        return true;
    }

    x = xOld;
    return false;
}



struct SAreturn{
    bool accepted = false;
    double L = 0.;
    double error = 0.;
};

template<typename Data, typename Func1, typename Func2>
SAreturn simulatedAnnealingStep(double temp, Data& par, Func1 L, Func2 T, Random& random){

    // can't do the following sice we need to save values of L (<H> in this case)
    //return metropolisStep(par, [&](Data par){return exp(- L(par) / temp);}, T, random);

    SAreturn ret;
    
    // we have to do this way since p (a lambda) in the Metropolis uses the istantaneus value of "par"
    // which changes during this function, and we want to compute <H> using the correct value:
    // "par" before applying T for <H> using old sigma and mu, and T(par) for the new sigma and mu
    vector<double> Lold = L(par);
    Data oldPar = par;
    par = T(par);
    vector<double> Lnew = L(par);

    ret.L = Lold[0];
    ret.error = Lold[1];

    double prob = exp(-(Lnew[0] - Lold[0]) / temp);

    if (random.Rannyu() <= min(1., prob)){
        ret.accepted = true;
        ret.L = Lnew[0];
        ret.error = Lnew[1];
    }
    else par = oldPar;

    return ret;
}

#endif