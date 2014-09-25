#ifndef BW_BOLTZMANN_WEIGHT_H_
#define BW_BOLTZMANN_WEIGHT_H_

#include "bw_types.h"
#include "bw_option_parser.h"
#include "bw_hdf5.h"

class Boltzmann_Weight : public Interact_H5
{
protected :
    bw_options options;
    std::vector<bw_datfile> list;
    std::vector<std::vector<double> > average;
    std::vector<std::vector<double> > stdev;
    bool xvg_numbered;
public :
    Boltzmann_Weight(const bw_options &option);
    void bw_read_filelist();
    void bw_parse_dat();
    void bw_read_dat(bw_datfile &file);
    void bw_write_dat(const bw_datfile file);
    void bw_calc_prob();
    void calc_average(h5_dat &prob, bw_datfile *list, const int &nexp, bool doHist);
    void bw_hist(const std::vector<std::vector<double> > &dat, int &index);
    ~Boltzmann_Weight(){};
};

#endif
