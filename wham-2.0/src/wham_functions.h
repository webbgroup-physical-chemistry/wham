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
    std::vector<std::vector<float> > trajectory;
    std::vector<float> c_major_omega, c_major_omega_transpose;
    std::vector<float> opt_trajectory, potential, probability;
public:
    void wham_init(t_wham args, t_options options);
    void TransposeOmega();
    std::vector<float> WhamStep(int point);
    float square_diff(std::vector<float> a, std::vector<float> b);
    void wham_pmf();
    void wham_prob();
    std::vector<float> PMF();
    std::vector<float> PROB();
};


#endif
