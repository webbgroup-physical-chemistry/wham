#ifndef WHAM_TRIG_H_
#define WHAM_TRIG_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>

#ifndef M_PI
#define M_PI = atan(1.)*4.
#endif

#ifndef RAD2DEG
#define RAD2DEG 180./M_PI
#endif

#ifndef DEG2RAD
#define DEG2RAD M_PI/180.
#endif


std::vector<float> angle2vector( float angle );

float periodic( float angle );

float delta_angle( float a,
                   float b);


template<class T> T vec_sum(std::vector<T> a)
{
    T result = 0;
    int asize = a.size();
    for (int i=0; i<asize; i++)
    {
        result += a[i];
    }
    return result;
}

template<class U, class R> U vec_dot(std::vector<U> a, std::vector<R> b)
{
    U product = 0;
    if (a.size() != b.size())
    {
        std::cerr << "\nERROR! Vector sizes do to match; cannot do dot product\n";
        std::exit(1);
    }
    int asize = a.size();
    for (int i=0; i<asize ; i++)
    {
        product += a[i]*b[i];
    }
    return product;
}

template<class V> void vec_normalize(std::vector<V> &a)
{
    V norm = vec_sum(a);
    int asize = a.size();
    for (int i=0; i<asize; i++)
    {
        a[i] /= norm;
    }
    return;
}


#endif