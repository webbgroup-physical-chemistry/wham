#include "bw_boltzmann_weight.h"
#include "bw_hdf5.h"

int main(int argc, char * argv[])
{
    bw_options opt;
    bw_option_parser(argc, argv, opt);
    
    Boltzmann_Weight results(opt);
//    results.bw_init(opt);
    return 0;
}
