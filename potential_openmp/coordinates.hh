#ifndef COORDINATES_HH
#define COORDINATES_HH

/* -------------------------------------------------------------------------- */
#include<stdio.h>
/* -------------------------------------------------------------------------- */
using namespace std;


struct coordinates {
double x, y, z;

coordinates(double u=0., double v=0., double w=0.);

coordinates operator- (coordinates u);

};

#endif /* COORDINATES_HH */

