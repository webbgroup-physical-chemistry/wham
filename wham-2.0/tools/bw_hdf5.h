#ifndef BW_HDF5_H_
#define BW_HDF5_H_

#include "bw_types.h"
#include "H5Cpp.h"
#include <sstream>

class Interact_H5
{
protected :
    bw_options options;
    std::string filename;
    H5::H5File *file;
public :
    Interact_H5(const bw_options &option);
    int h5_get_dataset(const std::string &path, h5_dat &data);
    void h5_bin_assignments(std::vector<bw_datfile> &list);
    // The write functions are all essentially the same and could/should
    // be consolidated at some point...
    void h5_write_dat(const bw_datfile &filedat);
    void h5_write_prob(const std::vector<h5_dat> &probs);
    void h5_write_hist(const std::vector<std::vector<double> > &dat, int &index);
    ~Interact_H5(){};
};

#endif
