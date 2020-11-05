/* -------------------------------------------------------------------------- */
#include "coordinates.hh"
/* -------------------------------------------------------------------------- */

coordinates::coordinates (double u, double v, double w)
{x = u;
y=v;
z=w;
}

coordinates coordinates::operator- (coordinates in) {
coordinates out;
out.x = x-in.x;
out.y = y-in.y;
out.z = z-in.z;
return out;
}

