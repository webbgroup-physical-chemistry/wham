#include "wham_trig.h"

std::vector<float> angle2vector( const float &angle )
{
    std::vector<float> v(2);
    v[0] = cos( angle * DEG2RAD );
    v[1] = sin( angle * DEG2RAD );
    return v;
}

float periodic( const float &angle )
{
    std::vector<float> v = angle2vector(angle);
    return atan2(v[1],v[0]) * RAD2DEG;
}

float delta_angle( const float &a, const float &b)
{
    std::vector<float> a1(2), b1(2);
    a1 = angle2vector(a);
    b1 = angle2vector(b);
    return acos(a1[0]*b1[0] + a1[1]*b1[1]) * RAD2DEG;
}



