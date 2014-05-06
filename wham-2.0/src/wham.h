#ifndef WHAM_H_
#define WHAM_H_
#include "parse_options.h"
#include "wham_types.h"
#include "wham_fileutils.h"
#include "wham_befeus.h"
#include "wham_functions.h"
#include <algorithm>
#include "wham_hdf5.h"

int cpp_wham(int argc, char *argv[]);

class WHAM
{
private:
    std::vector<TorsionExperiment> experiment;
    t_options options;
    int nstates, first_experiment, nexp, nzeros, nframes, last_frame;
    float invnstates;
    //std::vector<std::vector<TorsionExperiment> > group_traj;
    std::vector<t_experiment> group_traj;
    std::vector<t_map> map;
    std::vector<float> opt_prob, opt_pmf;
    t_wham wham_args;
    DOWHAM opttrajectory;
    Write_HDF5 outhdf5;
public:
    void cpp_wham_init(t_options option);
    void cpp_wham_read_experiments();
    void cpp_wham_find_nstates();
    void cpp_wham_initialize_vectors();
    void cpp_wham_group_ndim();
    void cpp_wham_assign_bins();
    void cpp_wham_count_bins();
    void cpp_wham_make_ND_map();
    void cpp_wham_mapOmega();
    void cpp_wham_doWHAM();
    void cpp_wham_conv();
    int map1d(std::vector<int> bins);
    
    //void cpp_wham_hdf5();


    int NFrames();
};



#endif