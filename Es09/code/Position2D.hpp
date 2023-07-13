#ifndef __Position2D__
#define __Position2D__

class Position2D{

public:
    double x, y;

    Position2D(){x = 0., y = 0.;}
    Position2D(double xx, double yy){x = xx; y = yy;}
    ~Position2D(){}

    Position2D operator+(const Position2D& other) const {return Position2D(x + other.x, y + other.y);}
    Position2D operator-(const Position2D& other) const {return Position2D(x - other.x, y - other.y);}
    double operator*(const Position2D& other) const {return x * other.x + y * other.y;}
};

#endif