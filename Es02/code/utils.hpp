#ifndef __utils__
#define __utils__

#include "Position.hpp"
#include "Random/random.h"

#include <cmath>

using namespace std;

double integrand(double x){
    return M_PI * 0.5 * cos(M_PI * 0.5 * x);
}

double calculateError(double sum, double sum2, int N){
    return (N > 1) ? sqrt(((sum2 - sum * sum / N) / N) / (N - 1)) : 0;
}

double probabilityDistribution(Random& random){
    return 1 + sqrt(1 - random.Rannyu());
}

// determines the axis of movement,
// then moves along it by a in the positive or negative direction
Position moveCubicLattice(double a, Random& random){
    Position p;
    int ax = random.Integer(0, 2);

    if (ax == 0){
        (random.Rannyu() < 0.5) ? p.x = a : p.x = -a;
    }
    else if (ax == 1){
        (random.Rannyu() < 0.5) ? p.y = a : p.y = -a;
    }
    else{
        (random.Rannyu() < 0.5) ? p.z = a : p.z = -a;
    }

    return p;
}

// move by a in a random direction in continuum 3D space
Position moveContinuum(double a, Random& random){
    Position p;

    // we generate a point on the surface of a unit sphere following what is explained
    // on the page https://mathworld.wolfram.com/SpherePointPicking.html (equations 9-11)
    // then, by multiplying it by a, we obtain the translation we need to apply to our point
    double x, y, r2;
    do{
        x = random.Rannyu(-1., 1.);
        y = random.Rannyu(-1., 1.);
        r2 = x * x + y * y;
    } while (r2 > 1.);

    p.x = a * 2. * x * sqrt(1. - x * x - y * y);
    p.y = a * 2. * y * sqrt(1. - x * x - y * y);
    p.z = a * (1. - 2. * (x * x + y * y));

    return p;
}

double randomSine(Random& random){
    double x, y, r2;
    do{
        x = random.Rannyu(-1., 1.);
        y = random.Rannyu(-1., 1.);
        r2 = x * x + y * y;
    } while (r2 > 1.);

    return y / sqrt(r2);
}

#endif