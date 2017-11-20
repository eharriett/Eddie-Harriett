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
