#ifndef BW_BOLTZMANN_WEIGHT_H_
#define BW_BOLTZMANN_WEIGHT_H_

#include "bw_types.h"
#include "bw_option_parser.h"
#include "bw_hdf5.h"

class Boltzmann_Weight
{
private :
    bw_options options;
    std::vector<bw_datfile> list;
    std::vector<std::vector<float> > average;
    std::vector<std::vector<float> > stdev;
    Interact_H5 h5file;
    bool xvg_numbered;
public :
    void bw_init(bw_options option);
    void bw_read_filelist();
    void bw_parse_dat();
    void bw_read_dat(bw_datfile &file);
    void bw_write_dat(bw_datfile file);
    void bw_calc_prob();
    void calc_average(h5_dat &prob, bw_datfile *list, int nexp);
};

#endif
