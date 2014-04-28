#include "wham.h"

int main(int argc, char * argv[])
{
    t_options opt;
    wham_options(argc, argv, opt);
    opt.ntraj = get_n_lines(opt.filelist);
    
    WHAM workflow;
    workflow.cpp_wham_init(opt);
    workflow.cpp_wham_conv();
    
    return 0;
}
