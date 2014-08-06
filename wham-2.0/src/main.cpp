#include "wham.h"

int main(int argc, char * argv[])
{
    t_options opt;
    wham_options(argc, argv, opt);
    opt.ntraj = get_n_lines(opt.filelist);
    
    WHAM workflow;
    workflow.cpp_wham_init(opt);
    workflow.cpp_wham_conv();
    /* Does not work currently.  Also not convinced it will save save
    if (opt.doConv)
    {
        cpp_wham_conv(workflow,opt);
    }
    */
    
    return 0;
}
