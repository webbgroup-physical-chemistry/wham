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
    double phi0, deltaphi, fc, t, beta, begin, end, bin_size;
    std::string xvgfile;
    std::vector<double> potential, omega;
    std::vector<double> coordinate;
    std::vector<int> bincounts;
    std::vector<t_bin> frame_bins;
    t_options options;
    bool precision_error;
public:
    void read_filelist_and_setup( std::string line,
                       t_options parser_options,
                       int previous_experiment,
                       int previous_dof );
    void read_trajectory();
    void build_histogram();
    void set_nbins();
    double trapezoidal_integration( double a, double b, int steps = 1000);
    double simpson_integration( double a, double b, int steps = 1000);
    double weight_dihedral_restraint_potential(double phi);
    void weight_dihedral_restraint_potentials();
    void assign_omega();
    void assign_NDbin(int frame, int bin);
    /* Return Functions for access class data */
    int Experiment();
    int DoF();
    int nBins();
    double T();
    double Phi0();
    double DPhi();
    double FC();
    std::vector<double> Potential();
    std::vector<double> Omega();
    std::vector<double> Coordinate();
    std::vector<int> BinCounts();
    std::vector<t_bin> FrameBins();
};


#endif