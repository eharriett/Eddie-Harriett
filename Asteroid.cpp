#include <iostream>
#include "Asteroid.h"
using namespace std;

const double G = .00006674;
const double t = .1;


Asteroid::Asteroid(double massin, double xin, double yin){
    public x = xin;
    public y = yin;
    public mass = massin;
}

void Asteroid::setForce(double forcexin, double forceyin){
    fx += forcexin;
    fy += forceyin;
}

double Asteroid::updatePosition(){
        xaccel = fx / mass;
        yaccel = fy / mass;
        xvel = xvel + xaccel * t;
        yvel = yvel + yaccel * t;
        x = x + xvel * t;
        y = y + yvel * t;
}

double Asteroid::dist(Asteroid a, Asteroid b)
    {
        return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
    }

double Asteroid::updateForce(Asteroid a, Asteroid b)
    {
        int slope = (a.y - b.y)/(a.x - b.x);
        if(slope > 1 || slope < -1){
            slope = slope - Trunc(slope);
        }
        int angle = atan(slope); 
        int forcex = (G * a.mass * b.mass)/(pow(dist(a,b),2)) * cos(angle);
        int forcey = (G * a.mass * b.mass)/(pow(dist(a,b),2)) * sin(angle);
        
        a.setForce(forcex, forcey);
        b.setForce(-forcex, forcey);
    }
