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
    std::vector<double> opt_prob, opt_pmf;
    std::string filename;
    H5::H5File *file;
public:
    void init(t_options option,
              t_wham args,
              int n_exp,
              int n_states,
              std::vector<double> prob,
              std::vector<double> pmf );
    void write_trj(std::vector<t_experiment> group_traj);
    void write_wham(std::vector<double> prob,
                    std::vector<double> pmf,
                    std::vector<int> counts,
                    std::vector<t_map> map,
                    std::string group,
                    int f0, int fN);
};