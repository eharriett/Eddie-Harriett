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
        Asteroid();
        Asteroid(double massin, double xin, double yin);
      //  ~Asteroid();
        void updatePosition();
        double massAst;
        double fx, fy;
        double x, y;
        double xvel, yvel, xaccel, yaccel;
        bool destroyed;
};

//default constructor
Asteroid::Asteroid(void){}

//default destructor
Asteroid::~Asteroid(void){}

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

void updatePlanet(Asteroid* a, Asteroid* b)
{
    double slope = (a->y - b->y)/(a->x - b->x);
    if(slope > 1 || slope < -1)
        slope = slope - trunc(slope);

    double angle = atan(slope); 
    double forcex = (G * a->massAst * b->massAst)/(pow(dist(a,b),2)) * cos(angle);
    double forcey = (G * a->massAst * b->massAst)/(pow(dist(a,b),2)) * sin(angle);

    a->fx = a->fx + forcex;
    a->fy = a->fy + forcey;
}

int main(int argc, char* argv[]) 
{
    if (argc < 6)
    {
        cout << "Not enough parameters" << endl;
        return 1;
    }

    for (int i = 0; i < argc; i++)
    {
        if (atoi(argv[i]) < 0)
        {
            cout << "Invalid input" <<endl;
            return 1;
        }
    }

    num_asteroids = atoi(argv[1]);
    num_iterations = atoi(argv[2]);
    num_planets = atoi(argv[3]);
    pos_ray = atoi(argv[4]);
    seed = atoi(argv[5]);
    
    vector<Asteroid>Asteroids(num_asteroids);
    vector<Asteroid>Planets(num_planets);
    
 //   Asteroid Asteroids[num_asteroids];
 //   Asteroid Planets[num_planets];

    default_random_engine re{seed}; 
    uniform_real_distribution<double> xdist{0.0, std::nextafter(spacewidth, std :: numeric_limits<double>::max())}; 
    uniform_real_distribution<double> ydist{0.0, std::nextafter( spaceheight, std :: numeric_limits<double>::max())}; 
    normal_distribution<double> mdist{mass, sdm};
    
    for (int i = 0; i < num_asteroids; i++)
    {
        Asteroid a(mdist(re), xdist(re), ydist(re));
        Asteroids[i] = a;
    }

    for (int i = 0; i < num_planets; i++)
    {
        if(i%4 == 0)
        {
            Asteroid a(mdist(re)*10, 0, ydist(re));
            Planets[i] = a;
        }
        else if(i%4 == 1)
        {
            Asteroid a(mdist(re)*10, xdist(re), 200);
            Planets[i] = a;
        }
        else if(i%4 == 2)
        {
            Asteroid a(mdist(re)*10, 200, ydist(re));
            Planets[i] = a;
        }
        else
        {
            Asteroid a(mdist(re)*10, xdist(re), 0);
            Planets[i] = a;
        }
    }

    for(int g = 0; g < num_iterations; g++)
    {
        for (int i = 0; i < num_asteroids-1; i++)
        {
            Asteroid* a;
            Asteroid* b;
            Asteroid* c;
            if(Asteroids[i].destroyed == false)
            {
                a = &Asteroids[i];
                for (int j = i+1; j < num_asteroids; j++)
                {
                    if(Asteroids[j].destroyed == false)
                    {
                        b = &Asteroids[j];
                        updateForce(a, b);
                    }
                }

                for (int h = 0; h < num_planets; h++)
                {
                    c = &Planets[h];
                    updatePlanet(a, c);
                }
            }
        }
        
        for(int i = 0; i < num_asteroids; i++)
        {
            if(Asteroids[i].destroyed == false)
            {
                updatePosition(&Asteroids[i]);

                if (Asteroids[i].x <= 0)
                    Asteroids[i].x = 2;

                if (Asteroids[i].y <= 0)
                    Asteroids[i].y = 2;

                if (Asteroids[i].x >= spacewidth)
                    Asteroids[i].x = spacewidth - 2;

                if (Asteroids[i].y >= spaceheight)
                    Asteroids[i].y = spaceheight - 2;

                //laser
                if (Asteroids[i].y >= pos_ray-2 && Asteroids[i].y <= pos_ray+2)
                        Asteroids[i].destroyed = true; 
            }
        }
    }
    
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
