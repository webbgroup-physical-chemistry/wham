#include "wham_trig.h"

std::vector<double> angle2vector( double angle )
{
    std::vector<double> v(2);
    v[0] = cos( angle * DEG2RAD );
    v[1] = sin( angle * DEG2RAD );
    return v;
}

double periodic( double angle )
{
    std::vector<double> v = angle2vector(angle);
    return atan2(v[1],v[0]) * RAD2DEG;
}

double delta_angle( double a, double b)
{
    std::vector<double> a1(2), b1(2);
    a1 = angle2vector(a);
    b1 = angle2vector(b);
    return acos(a1[0]*b1[0] + a1[1]*b1[1]) * RAD2DEG;
}



