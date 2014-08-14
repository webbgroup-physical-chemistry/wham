#ifndef BW_HDF5_H_
#define BW_HDF5_H_

#include "bw_types.h"
#include "H5Cpp.h"
#include <sstream>

class Interact_H5
{
private :
    bw_options options;
    std::string filename;
    H5::H5File *file;
public :
    void h5_init(const bw_options &option);
    void h5_write_dat(const bw_datfile &filedat);
    int h5_get_dataset(const std::string &path, h5_dat &data);
    void h5_bin_assignments(std::vector<bw_datfile> &list);
    void h5_write_prob(const std::vector<h5_dat> &probs);
};

#endif
