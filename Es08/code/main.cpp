#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "utils.hpp"
#include "Random/random.h"

using namespace std;

int main(){
    double const hbar = 1., m = 1.;
    double stepMetro, stepSA, sigma0, mu0, temp0, dtemp;
    int nBlocks, nSteps, nSpT; // steps per temperature
    int accepted = 0, attempted = 0;
    int acceptedSA = 0, attemptedSA = 0;

    Random random;

    ifstream in("input.in");
    in >> sigma0 >> mu0;
    in >> stepMetro >> nBlocks >> nSteps;
    in >> stepSA >> temp0 >> dtemp >> nSpT;
    in.close();

    vector<double> params = {sigma0, mu0};

    double x = generateStartingPoint(sigma0, mu0, random);

    // transition probability for computing <H>
    auto T = [&stepMetro, &random](double x){
        return x + stepMetro * (random.Rannyu() - 0.5);
    };

    // in the Metropolis algorithm we are only interested in ratios of probabilities
    // so we can ignore the normalization
    auto p = [&params](double x){
        double mod = psiT(x, params[0], params[1]);
        return mod * mod;
    };

    // transition probability for the SA algorithm (for mu and sigma)
    auto Tvec = [&stepSA, &random](vector<double> x){
        vector<double> xNew;
        for (int i = 0; i < (int) x.size(); i++) xNew.push_back(x[i] + stepSA * (random.Rannyu() - 0.5));
        return xNew;
    };

    // energy at x
    auto Eloc = [&hbar, &m](double x, vector<double> par){
        double sigma = par[0], mu = par[1];
        return (- hbar * hbar * 0.5 / m * d2psiT(x, sigma, mu)) / psiT(x, sigma, mu) + potential(x);
    };

    // computes <H> for one value of sigma and mu
    auto Hvar = [&](vector<double> par){
        double sumTot = 0., sumTot2 = 0.;

        for (int iBlock = 1; iBlock <= nBlocks; iBlock++){
            double sum = 0.;

            for (int iStep = 0; iStep < nSteps; iStep++){

                if (metropolisStep(x, p, T, random)) accepted++;
                attempted++;
                
                sum += Eloc(x, par);
            }

            double s = sum / nSteps;
            sumTot += s;
            sumTot2 += s * s;
        }

        return vector<double>{sumTot / nBlocks, error(sumTot, sumTot2, nBlocks)};
    };



    // temperature update laws
    auto updateTemperature = [&](double T){
        return T - dtemp;
    };

    auto stepsPerTemperature = [&](double T){
        return (int) nSpT;
    };



    // update sigma and mu using the simulated annealing algorithm
    const int wd = 12;

    cout << "starting parameters: temperature = " << temp0 << "\tsigma = " << sigma0 << "\tmu = " << mu0 << endl << endl;
    cout << setw(wd) << "temperature" << setw(wd) << "sigma" << setw(wd) << "mu" << setw(wd) << "<H>" << setw(wd) << "metro rate" << setw(wd) << "SA rate" << endl;

    ofstream outSM("sigma_mu.out");
    ofstream outH("H.out");

    for (double temp = temp0; temp > 0 + dtemp * 0.1; temp = updateTemperature(temp)){
        
        accepted = 0, attempted = 0;
        acceptedSA = 0, attemptedSA = 0;
        SAreturn SA;

        for (int i = 0; i < stepsPerTemperature(temp); i++){

            SA = simulatedAnnealingStep(temp, params, Hvar, Tvec, random);
            outH << temp << "\t" << SA.L << "\t" << SA.error << endl;
            if (SA.accepted) acceptedSA++;
            attemptedSA++;
        }

        cout << setw(wd) << temp << setw(wd) << params[0] << setw(wd) << params[1] << setw(wd) << SA.L << setw(wd) << (double) accepted / attempted << setw(wd) << (double) acceptedSA / attemptedSA << endl;
        outSM << temp << "\t" << params[0] << "\t" << params[1] << endl;
    }

    outSM.close();
    outH.close();

    cout << "annealing done" << endl;



    cout << "computing final value of <H>" << endl << endl;

    // estimate <H> with current values of sigma and mu using the metropolis algorithm
    ofstream out("hamiltonian.out");
    ofstream points("points.out");

    double sum = 0., sumOverBlocks = 0., sum2OverBlocks = 0.;

    for (int iBlock = 1; iBlock <= nBlocks; iBlock++){
        sum = 0.;
        accepted = 0, attempted = 0;

        for (int iStep = 0; iStep < nSteps; iStep++){
            if (metropolisStep(x, p, T, random)) accepted++;
            attempted++;
            
            sum += Eloc(x, params);

            points << x << endl;
        }

        double s = sum / nSteps;
        sumOverBlocks += s;
        sum2OverBlocks += s * s;

        out << s << "\t" << sumOverBlocks / iBlock << "\t" << error(sumOverBlocks, sum2OverBlocks, iBlock) << endl;

        cout << "Block " << iBlock << " done\n";
        cout << "acceptance rate: " << (double)accepted/(double)attempted << endl << endl;
    }

    out.close();
    points.close();



    return 0;
}