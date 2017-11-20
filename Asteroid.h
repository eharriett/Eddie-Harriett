#ifndef _ASTEROID_H_
#define _ASTEROID_H

class Asteroid
{
public: 
	Asteroid(double massin, double xin, double yin);
	void setForce(double fxin, double fyin);
	double updatePosition();
	double mass;
	double fx = 0, fy = 0;
	double x, y, xvel = 0, yvel = 0, xaccel = 0, yaccel = 0;

};

#endif
