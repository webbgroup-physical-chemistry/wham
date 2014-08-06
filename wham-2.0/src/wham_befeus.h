#ifndef WHAM_BEFEUS_H_
#define WHAM_BEFEUS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include "wham_trig.h"
#include "wham_types.h"

class TorsionExperiment
{
private:
    int experiment, dof, nbins;
    float phi0, deltaphi, fc, t, beta, begin, end, bin_size;
    std::string xvgfile;
    std::vector<float> potential;
    std::vector<double> omega;
    std::vector<float> coordinate;
    std::vector<int> bincounts;
    std::vector<t_bin> frame_bins;
    t_options options;
    bool precision_error;
public:
    void read_filelist( std::string line,
                        t_options &parser_options,
                        int previous_experiment,
                        int previous_dof );
    void setup(t_options &parser_options);
    void read_trajectory();
    void build_histogram();
    void set_nbins();
    double trapezoidal_integration( float a, float b, int steps = 1000);
    double simpson_integration( float a, float b, int steps = 1000);
    float weight_dihedral_restraint_potential(float phi);
    void weight_dihedral_restraint_potentials();
    void assign_omega();
    void assign_NDbin(int frame, int bin);
    /* Return Functions for access class data */
    int Experiment();
    int DoF();
    int nBins();
    float T();
    float Phi0();
    float DPhi();
    float FC();
    std::vector<float> Potential();
    std::vector<double> Omega();
    std::vector<float> Coordinate();
    std::vector<int> BinCounts();
    std::vector<t_bin> FrameBins();
};


#endif
