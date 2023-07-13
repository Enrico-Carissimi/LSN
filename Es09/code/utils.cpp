#include <cmath>
#include "utils.hpp"

using namespace std;



vector<double> generateRandom(int l, double min, double max, Random& random){
    if (min > max) swap<double>(min, max);

    vector<double> v;

    for (int i = 0; i < l; i++){
        v.push_back(random.Rannyu(min, max));
    }

    return v;
}



vector<Position2D> generateRandom2D(int l, double r, Random& random){
    vector<Position2D> v;

    for (int i = 0; i < l; i++){
        double angle = random.Rannyu(0., 2 * M_PI);
        v.push_back(Position2D(r * cos(angle), r * sin(angle)));
    }

    return v;
}



vector<Position2D> generateRandom2D(int l, double xmin, double xmax, double ymin, double ymax, Random& random){
    vector<double> x = generateRandom(l, xmin, xmax, random);
    vector<double> y = generateRandom(l, ymin, ymax, random);

    vector<Position2D> v;

    for (int i = 0; i < l; i++){
        v.push_back(Position2D(x[i], y[i]));
    }

    return v;
}



// generates cartesian coordinates using angles, the radius can be scaled after
vector<Position2D> generateOnCircle(vector<double> angles){
    vector<Position2D> v;

    for (int i = 0; i < (int)angles.size(); i++){
        v.push_back(Position2D(cos(angles[i]), sin(angles[i])));
    }

    return v;
}



double distanceOnCircle(double r, double a, double b){
    double d = abs(a - b);
    if (d > M_PI) d = 2. * M_PI - d;

    return r * d;
}



// the following two functions can be optimized more: we compute N^2 distances
// but with N elements, there are actually N(N-1)/2 couples
// the distance between element i and j is stored in the element i*size+j
vector<double> computeDistancesCircle(int size, double r, vector<double> angles){
    vector<double> distances;

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            distances.push_back(distanceOnCircle(r, angles[i], angles[j]));
        }
    }

    return distances;
}



vector<double> computeDistancesSquare(int size, vector<Position2D> pos){
    vector<double> distances;

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            distances.push_back(distance2(pos[i], pos[j]));
        }
    }

    return distances;
}



void printPositions(string fileName, vector<Position2D> v){
    ofstream out(fileName);

    for (int i = 0; i < (int)v.size(); i++){
        out << v[i].x << "\t" << v[i].y << endl;
    }

    out.close();
}