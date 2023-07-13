#include <cmath>
#include "utils.hpp"

using namespace std;



vector<Position2D> readCapitals(){
    vector<Position2D> coords;

    string buffer;
    double longitude, latitude;
    ifstream File("American_capitals.dat");

    // skip first line
    File >> buffer >> buffer >> buffer >> buffer;

    // read content (only coordinates)
    while (File >> buffer >> buffer >> longitude >> latitude){
        coords.push_back(Position2D(latitude, longitude));
    }

    // convert latitude and longitude in radians
    for (int i = 0; i < (int)coords.size(); i++){
        coords[i].x *= M_PI / 180.;
        coords[i].y *= M_PI / 180.;
    }

    return coords;
}



// we compute distances as the shortest path on the sphere
double distanceOnSphere(double r, Position2D a, Position2D b){

    // x is lat, y is long
    double d = sin(a.x) * sin(b.x) + cos(a.x) * cos(b.x) * cos(a.y - b.y);
    
    return r * acos(d);
}



// the following function can be optimized more: we compute N^2 distances
// but with N elements, there are actually N(N-1)/2 couples
// the distance between element i and j is stored in the element i*size+j
vector<double> computeDistances(int size, double r, vector<Position2D> angles){
    vector<double> distances;

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            distances.push_back(distanceOnSphere(r, angles[i], angles[j]));
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