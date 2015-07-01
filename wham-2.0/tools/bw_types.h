#ifndef BW_TYPES_H_
#define BW_TYPES_H_

#include <vector>
#include <string>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <fstream>
#include <iostream>

#ifndef M_PI
#define M_PI = atan(1.)*4.
#endif

#ifndef RAD2DEG
#define RAD2DEG 180./M_PI
#endif

#ifndef DEG2RAD
#define DEG2RAD M_PI/180.
#endif

inline double MAX(double a, double b) {
    return (a>b)? a : b;
}
inline double MIN(double a, double b) {
    return (a<b)? a : b;
}

struct bw_options
{
    std::string xvglist;
    std::string hdf5file;
    std::vector<std::string> datnames;
    std::vector<std::string> datunits;
    int nitems;
    int frameStep;
    int histbins;
    int seed;
    bool doseed;
    bool isangle;
    bool bVerbose;
    double nbootstrap;
    bool bootstrap;
    bool doWrite;
};

struct bw_datfile
{
    int experiment;
    std::string filename;
    std::vector<std::vector<double> > dat;
    std::vector<int> frameN;
    std::vector<int> bin;
};

struct h5_dat
{
    std::vector<int> span;
    std::vector<double> bin_prob;
    std::vector<double> avg;
    std::vector<double> stdev;
};



#endif
