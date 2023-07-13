#ifndef __Position__
#define __Position___

class Position{
public:
    
    double x, y, z;
    
    Position(){x = 0., y = 0., z = 0.;}
    Position(double x, double y, double z){setPosition(x, y, z);}
    Position(const Position& p){setPosition(p);}

    ~Position(){}
    
    void setPosition(double X, double Y, double Z){
        x = X, y = Y, z = Z;
    }
    void setPosition(const Position& p){
        x = p.x, y = p.y, z = p.z;
    }

    double distance2(double x2, double y2, double z2){
        return (x2 - x) * (x2 - x) + (y2 - y) * (y2 - y) + (z2 - z) * (z2 - z);
    }
    double distance2(const Position& p){
        return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z) * (p.z - z);
    }

    void move(double i, double j, double k){
        x += i, y += j, z += k;
    }
    //moving a position by another "position" doesn't really make sense
    //we should create another "vector3d" class but it's not needed in this case
    void move(const Position& p){
        x += p.x, y += p.y, z += p.z;
    }
};

#endif