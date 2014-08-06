#ifndef WHAM_TYPES_H_
#define WHAM_TYPES_H_

#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <limits>
#include <ctime>
#include <cstdlib>
#include <stdio.h>

#ifndef KB
#define KB 0.00813144621 // kJ/(mol*K)
#endif

struct t_options
{
    std::string filelist;
    std::string outname;
    int f0;
    int fN;
    float t;
    std::vector<float> x0;
    std::vector<float> xN;
    std::vector<float> b;
    int ndof;
    float tol;
    float iter;
    int ntraj;
    bool bVerbose;
    bool doConv;
    int convStep;
    int nexceed;
};

struct t_bin
{
    int bin1D;
    int binND;
};

struct t_map
{
    std::vector<int> bin1D;
};

struct t_wham
{
    int nstates;
    std::vector<int> sample;
    std::vector<int> counts;
    std::vector<float> t;
    std::vector<std::vector<double> > omegas;
};

struct t_experiment
{
    int experiment;
    float t;
    // [dofM]
    std::vector<int> nbins;
    std::vector<float> phi0, deltaphi, fc;
    // [frameN][dofM]
    std::vector<std::vector<float> > potential, coordinate;
    std::vector<std::vector<t_bin> > frame_bin;
    // [binL][dofM]
    std::vector<std::vector<double> > omega;
    std::vector<std::vector<int> > bincounts;
};

#endif
