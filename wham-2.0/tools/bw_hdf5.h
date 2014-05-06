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
    void h5_init(bw_options option);
    void h5_write_dat(bw_datfile filedat);
    int h5_get_dataset(std::string path, h5_dat &data);
    void h5_bin_assignments(std::vector<bw_datfile> &list);
    void h5_write_prob(std::vector<h5_dat> probs);
};

#endif
