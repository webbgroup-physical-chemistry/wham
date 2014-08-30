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
#ifdef _OPENMP
#include <omp.h>
#endif

class DOWHAM
{
private:
    t_wham wham_args;
    t_options wham_options;
    int nexperiments;
    //std::vector<std::vector<float> > trajectory;
    std::vector<float> c_major_omega, c_major_omega_transpose;
    std::vector<float> opt_trajectory, potential, probability;
public:
    void wham_init(const t_wham &args, const t_options &options);
    void TransposeOmega();
    void WhamStep(std::vector<float> &step);
    float square_diff(const float *a, const float *b, const int &sizea, const int &sizeb);
    void wham_pmf();
    void wham_prob();
    std::vector<float> PMF();
    std::vector<float> PROB();
};


#endif
