#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <fstream>
#include <iomanip>

using namespace std;

//declaration of constants
int num_asteroids, num_iterations, num_planets, seed;
const double G = .00006674;
const int spacewidth = 200;
const int spaceheight = 200;
const double mass = 1000;
const int sdm = 50;
const double t = .1;
const double ray_width = 2;
double pos_ray;

//class for asteroids and planets
class Asteroid{
    public: 
        Asteroid(double massin, double xin, double yin);
        void updatePosition();
        double massAst;
        double fx, fy;
        double x, y;
        double xvel, yvel, xaccel, yaccel;
        bool destroyed;
};

//constructor
Asteroid::Asteroid(double massin, double xin, double yin)
{
    x = xin;
    y = yin;
    massAst = massin;
    fx = 0; fy = 0; xvel = 0; yvel = 0; xaccel = 0; yaccel = 0;
    destroyed = false;
}

//updates the all values for the asteroids
void updatePosition(Asteroid* a)
{
    if(a->fx > 200){
        a->fx = 200;
    }
    if(a->fy > 200){
        a->fy = 200;
    }
    
    //updates acceleration
    a->xaccel = a->xaccel + a->fx / a->massAst; //a = f/m
    a->yaccel = a->yaccel + a->fy / a->massAst;
    
    //updates velocity
    a->xvel = a->xvel + a->xaccel * t; //v = at
    a->yvel = a->yvel + a->yaccel * t;
    
    //updates position
    a->x = a->x + a->xvel * t; //p = vt
    a->y = a->y + a->yvel * t;

    //resets the forces for the next iteration
    a->fx = 0;
    a->fy = 0;
}

//calculates the distance betweek 2 objects
double dist(Asteroid *a, Asteroid *b)
{
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
}

//calculates the force between two objects and updates both
void updateForce(Asteroid* a, Asteroid* b)
{
    //finds slope
    double slope = (a->y - b->y)/(a->x - b->x);
    if(slope > 1 || slope < -1)
        slope = slope - trunc(slope);
    
    //finds angle 
    double angle = atan(slope); 
    double forcex = (G * a->massAst * b->massAst)/(pow(dist(a,b),2)) * cos(angle);
    double forcey = (G * a->massAst * b->massAst)/(pow(dist(a,b),2)) * sin(angle);

    //updates forces for object a
    a->fx = a->fx + forcex;
    a->fy = a->fy + forcey;

    //updates forces for object b
    b->fx = b->fx - forcex;
    b->fy = b->fy - forcey;
}

//updates asteroid values based on the planets
void updatePlanet(Asteroid* a, Asteroid* b)
{
    //finds slope
    double slope = (a->y - b->y)/(a->x - b->x);
    if(slope > 1 || slope < -1)
        slope = slope - trunc(slope);

    //finds angle
    double angle = atan(slope); 
    double forcex = (G * a->massAst * b->massAst)/(pow(dist(a,b),2)) * cos(angle);
    double forcey = (G * a->massAst * b->massAst)/(pow(dist(a,b),2)) * sin(angle);

    //updates forces for the asteroid
    a->fx = a->fx + forcex;
    a->fy = a->fy + forcey;
}

int main(int argc, char* argv[]) 
{
    //program isn't executed with all parameters
    if (argc < 6)
    {
        cout << "Not enough parameters" << endl;
        return 1;
    }

    //checks if any inputs are less than 0
    for (int i = 0; i < argc; i++)
    {
        if (atoi(argv[i]) < 0)
        {
            cout << "Invalid input" <<endl;
            return 1;
        }
    }

    //sets inputted values to corresponding variables
    num_asteroids = atoi(argv[1]);
    num_iterations = atoi(argv[2]);
    num_planets = atoi(argv[3]);
    pos_ray = atoi(argv[4]);
    seed = atoi(argv[5]);
    
    //creates vectors for asteroids and planets
    vector<Asteroid>Asteroids;
    vector<Asteroid>Planets;
    
    //creates random numbers
    default_random_engine re{seed}; 
    uniform_real_distribution<double> xdist{0.0, std::nextafter(spacewidth, std :: numeric_limits<double>::max())}; 
    uniform_real_distribution<double> ydist{0.0, std::nextafter( spaceheight, std :: numeric_limits<double>::max())}; 
    normal_distribution<double> mdist{mass, sdm};
    
    #pragma omp parallel num_threads(4)
    #pragma omp for
    //assigns random values for the asteroids
    for (int i = 0; i < num_asteroids; i++)
    {
        Asteroid a(mdist(re), xdist(re), ydist(re));
        Asteroids.push_back(a);
    }

    #pragma omp parallel num_threads(4)
    #pragma omp for    
    //assigns random values for the planets
    for (int i = 0; i < num_planets; i++)
    {
        if(i%4 == 0)
        {
            Asteroid a(mdist(re)*10, 0, ydist(re));
            Planets.push_back(a);
        }
        else if(i%4 == 1)
        {
            Asteroid a(mdist(re)*10, xdist(re), 200);
            Planets.push_back(a);
        }
        else if(i%4 == 2)
        {
            Asteroid a(mdist(re)*10, 200, ydist(re));
            Planets.push_back(a);
        }
        else
        {
            Asteroid a(mdist(re)*10, xdist(re), 0);
            Planets.push_back(a);
        }
    }

    //loop through iterations
    for(int g = 0; g < num_iterations; g++)
    {
        #pragma omp parallel num_threads(4)
        #pragma omp for
        //loop through all asteroids in order to update the force
        for (int i = 0; i < num_asteroids-1; i++)
        {
            //creates pointers to an asteroid so we can pass them to other functions
            Asteroid* a;
            Asteroid* b;
            Asteroid* c;
            
            //checks if asteroid a is still valid
            if (Asteroids[i].destroyed == false)
            {
                a = &Asteroids[i];
                
                //checks asteroid a against all other asteroids
                for (int j = i+1; j < num_asteroids; j++)
                {
                    //checks is asteroid b is still valid
                    if(Asteroids[j].destroyed == false)
                    {
                        b = &Asteroids[j];
                        updateForce(a, b);
                    }
                }

                //checks asteroid a against all planets
                for (int h = 0; h < num_planets; h++)
                {
                    c = &Planets[h];
                    updatePlanet(a, c);
                }
            }
        }
        
        #pragma omp parallel num_threads(4)
        #pragma omp for
        //loop through all asteroids in order to update the positions
        for(int i = 0; i < num_asteroids; i++)
        {
            if(Asteroids[i].destroyed == false)
            {
                updatePosition(&Asteroids[i]);

                //if the asteroid's x or y position is within 2 from the edge, set it to 2
                if (Asteroids[i].x < 2){
                    Asteroids[i].x = 2;
                    Asteroids[i].xvel = -Asteroids[i].xvel;
                }

                if (Asteroids[i].y < 2){
                    Asteroids[i].y = 2;
                    Asteroids[i].yvel = -Asteroids[i].yvel;
                }

                if (Asteroids[i].x > spacewidth - 2){
                    Asteroids[i].x = spacewidth - 2;
                    Asteroids[i].xvel = -Asteroids[i].xvel;
                }                    

                if (Asteroids[i].y > spaceheight - 2){
                    Asteroids[i].y = spaceheight - 2;
                    Asteroids[i].yvel = -Asteroids[i].yvel;
                }

                //checks if asteroid hit the laser
                if (Asteroids[i].y >= pos_ray-2 && Asteroids[i].y <= pos_ray+2)
                        Asteroids[i].destroyed = true; 
            }
        }
    }
    
    //writes to file called out.txt
    ofstream inFile;
    inFile.open("out.txt");
    if (inFile.is_open())
    {
        for (int g = 0; g < num_asteroids; g++)
        {
            if (Asteroids[g].destroyed == false)
            {
                cout << Asteroids[g].y << endl;
                inFile << fixed << setprecision(3) << Asteroids[g].y << " " << Asteroids[g].x << " " << Asteroids[g].xvel << " " << Asteroids[g].yvel << " " << Asteroids[g].massAst << "\n";
            }
        }
        inFile.close();
    }
    else cout << "Unable to open file";
    
    return 0;
}
