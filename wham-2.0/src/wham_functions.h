#ifndef WHAM_FUNCTIONS_H_
#define WHAM_FUNCTIONS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iterator>
#include <cstdlib>
#include "wham_types.h"
#include "wham_trig.h"

class DOWHAM
{
private:
    t_wham wham_args;
    t_options wham_options;
    int nexperiments;
    std::vector<std::vector<double> > trajectory, omegaT;
    std::vector<double> opt_trajectory, potential, probability;
public:
    void wham_init(t_wham args, t_options options);
    void TransposeOmega();
    std::vector<double> WhamStep(int point);
    double square_diff(std::vector<double> a, std::vector<double> b);
    void wham_pmf();
    void wham_prob();
    std::vector<double> PMF();
    std::vector<double> PROB();
};


#endif