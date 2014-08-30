#ifndef WHAM_HDF5_H_
#define WHAM_HDF5_H_

#include "wham_types.h"
#include "wham_befeus.h"
#include "H5Cpp.h"
#include <sstream>
#include <cmath>

class Write_HDF5
{
private:
    t_options options;
    t_wham wham_args;
    int nexp, nstates;
    std::vector<float> opt_prob, opt_pmf;
    std::string filename;
    H5::H5File *file;
public:
    void init(const t_options &option,
              const t_wham &args,
              const int &n_exp,
              const int &n_states,
              const std::vector<float> &prob,
              const std::vector<float> &pmf );
    void write_trj(const std::vector<t_experiment> &group_traj);
    void write_wham(const std::vector<float> &prob,
                    const std::vector<float> &pmf,
                    const std::vector<int> &counts,
                    const std::vector<t_map> &map,
                    const std::string &group,
                    const int &f0, const int &fN);
};

#endif
