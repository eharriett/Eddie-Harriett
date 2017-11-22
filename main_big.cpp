/*#include <nasteroids_seq.h>
#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include "Asteroid.h"

using namespace std;

int num_asteroids, num_iterations, num_planets, seed;
const int height = 2; 
const int width = 2;
const int mass = 1000;
const int sdm = 50;
double pos_ray;

int main() 
{
    cout << "Please enter the number asteroids to be simulated, greater or equal to 0: " << endl; 
    cout << "Please enter the number of iterations that will be simulated, greater or equal to 0: " << endl;
    cout << "Please enter the number of planets that will be simulated, greater or equal to 0: " << endl;
    cout << "Please enter the position of the ray in space: " << endl;
    cout << "Please enter a positive number for the random number generator: " << endl;
    cin >> num_asteroids >> num_iterations >> num_planets >> pos_ray >> seed;

    if (cin.fail() || num_asteroids < 0 || num_planets < 0 || num_iterations < 0 || seed < 0)
    {
    	cout << "Input Failure." << endl;
    	exit;
    }
    
    vector <Asteroid> Asteroids;
    vector <Asteroid> Planets;

    default_random_engine re{seed}; 
	uniform_real_distribution<double> xdist{0.0, std::nextafter( width, std :: numeric_limits<double>::max())}; 
	uniform_real_distribution<double> ydist{0.0, std::nextafter( height, std :: numeric_limits<double>::max())}; 
	normal_distribution<double> mdist{mass, sdm};
    
    for (int i = 0; i < (num_asteroids); i++)
    {
        Asteroids[i].x = xdist(re);
        Asteroids[i].y = ydist(re);
        Asteroids[i].mass = mdist(re);
    }

    for (int i = 0; i < (num_planets); i++)
    {
        Planets[i].x = xdist(re);
        Asteroids[i].y = ydist(re);
        Asteroids[i].mass = mdist(re) * 10;
    }

    for (int i = 0; i < (num_asteroids-1); i++)
    {
        for (int j = 0; i < num_asteroids; j++)
        {
            Asteroid a = Asteroids[i];
            Asteroid b = Asteroids[j];
            double dist = Asteroid.distance(a, b);

            // long and annoying parameters for space dimensions
            if (a.x <= 0)
            {
                a.x = 2;
            }

            if (a.y <= 0)
            {
                a.y = 2;
            }

            if (a.x >= Asteroid.width)
            {
                a.x = Asteroid.width - 2;
            }

            if (a.y >= Asteroid.width)
            {
                a.y Asteroid.width - 2;
            }

            if (dist > 2)
            {
                Asteroid.updateForce(a, b);
                // this is definitely fucked up
                a.updatePosition();

            }



        }
    }

    
    return 0;
}
*/

#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>

using namespace std;

int num_asteroids, num_iterations, num_planets, seed;
const double G = .00006674;
const int spacewidth = 200;
const int spaceheight = 200;
const double height = 2; 
const double width = 2;
const double mass = 1000;
const int sdm = 50;
const double t = .1;
double pos_ray;

class Asteroid{
    public: 
        Asteroid();
        Asteroid(double massin, double xin, double yin);
        void setForce(double fxin, double fyin);
        void updatePosition();
        double massAst;
        double fx, fy;
        double x, y;
        double xvel, yvel, xaccel, yaccel;
};
Asteroid::Asteroid(void){

}
Asteroid::Asteroid(double massin, double xin, double yin){
    x = xin;
    y = yin;
    massAst = massin;
   // fx = 0; fy = 0; xvel = 0; yvel = 0; xaccel = 0; yaccel = 0;
}

void setForce(Asteroid &a, double forcex, double forcey){
    cout << a.fx << " fx" << endl;
    a.fx = a.fx + forcex;
    cout << a.fx << " fx after" << endl;
    a.fy = a.fy + forcey;
}

void updatePosition(Asteroid &a){
    a.xaccel = a.fx / a.massAst;
    a.yaccel = a.fy / a.massAst;
    a.xvel = a.xvel + a.xaccel * t;
    a.yvel = a.yvel + a.yaccel * t;
    a.x = a.x + a.xvel * t;
    a.y = a.y + a.yvel * t;
}

double dist(Asteroid &a, Asteroid &b)
{
    double dista = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
    return dista;

}

double updateForce(Asteroid &a, Asteroid &b)
{
    double slope = (a.y - b.y)/(a.x - b.x);
    if(slope > 1 || slope < -1){
        slope = slope - trunc(slope);
    }
    double angle = atan(slope); 
    double forcex = (G * a.massAst * b.massAst)/(pow(dist(a,b),2)) * cos(angle);
    double forcey = (G * a.massAst * b.massAst)/(pow(dist(a,b),2)) * sin(angle);

    //setForce(*a, forcex, forcey);
    //setForce(*b, -forcex, -forcey);
    cout << a.fx << " fx" << endl;
    a.fx = a.fx + forcex;
    cout << a.fx << " fx after" << endl;
    a.fy = a.fy + forcey;

    b.fx = b.fx - forcex;
    b.fy = b.fy - forcey;
}

int main() 
{
    cout << "Please enter the number asteroids to be simulated, greater or equal to 0: " << endl; 
    cout << "Please enter the number of iterations that will be simulated, greater or equal to 0: " << endl;
    cout << "Please enter the number of planets that will be simulated, greater or equal to 0: " << endl;
    cout << "Please enter the position of the ray in space: " << endl;
    cout << "Please enter a positive number for the random number generator: " << endl;
    cin >> num_asteroids >> num_iterations >> num_planets >> pos_ray >> seed;

    if (cin.fail() || num_asteroids < 0 || num_planets < 0 || num_iterations < 0 || seed < 0)
    {
        cout << "Input Failure." << endl;
        //exit;
    }
    
    Asteroid Asteroids[num_asteroids];
    //Planet Planets[some number];

    default_random_engine re{seed}; 
    uniform_real_distribution<double> xdist{0.0, std::nextafter(spacewidth, std :: numeric_limits<double>::max())}; 
    uniform_real_distribution<double> ydist{0.0, std::nextafter( spaceheight, std :: numeric_limits<double>::max())}; 
    normal_distribution<double> mdist{mass, sdm};
    
    for (int i = 0; i < num_asteroids; i++)
    {
        Asteroid a(mdist(re), xdist(re), ydist(re));
        Asteroids[i] = a;
    }

   /* for (int i = 0; i < (num_planets); i++)
    {
        Planets[i].x = xdist(re);
        Asteroids[i].y = ydist(re);
        Asteroids[i].mass = mdist(re) * 10;
    } */
    for(int g = 0; g < num_iterations; g++){
        for (int i = 0; i < num_asteroids-1; i++)
        {
            for (int j = i+1; j < num_asteroids; j++)
            {
                Asteroid a = Asteroids[i];
                Asteroid b = Asteroids[j];

                updateForce(a, b);
            }
        }
        for(int i = 0; i < num_asteroids; i++){
            Asteroid a = Asteroids[i];
            updatePosition(a);
            if (a.x <= 0)
            {
                a.x = 2;
            }

            if (a.y <= 0)
            {
                a.y = 2;
            }

            if (a.x >= spacewidth)
            {
                a.x = spacewidth - 2;
            }

            if (a.y >= spaceheight)
            {
                a.y = spaceheight - 2;
            }
        }

        cout << Asteroids[1].fx << endl;

    }
    
    return 0;
}
