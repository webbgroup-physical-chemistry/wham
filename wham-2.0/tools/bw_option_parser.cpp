#include "bw_option_parser.h"

void bw_option_parser(int argc, char *argv[], bw_options &options)
{
    namespace po = boost::program_options;
    
    const char *description =
    {
        "\tUse to calculate the boltzmann weight average value based on the probability distribution generated from bin/wham.\n"
    };
    /*options.seed = 47;
    options.doseed = false;*/
    options.bootstrap = false;
    options.nbootstrap = 0;
    options.bVerbose = true;
    options.isangle = false;
    try
    {
        po::options_description desc(description);
        desc.add_options()
            ("help,h","")
            ("xvglist,f",po::value<std::string>(&options.xvglist)->required(),
                "File containing the list of column-data files for each trajectory.  Each column is a different property.  There must be an equal number of --datatype entries as columns.  Required")
            ("h5file,p",po::value<std::string>(&options.hdf5file)->required(),
                "HDF5 file for the system being examined.  Required")
            ("datatype,d",po::value<std::vector<std::string> >(&options.datnames)->required(),
                "Data type names for appending to HDF5 file.  Do not use spaces.  Required")
            ("dataunit,u",po::value<std::vector<std::string> >(&options.datunits)->required(),
                "Data type units for appending to HDF5 file.  Do not use spaces.  Required")
            /*("seed,s",po::value<int>(&options.seed),
                "Random seed to use for bootstrapping.  Default: 47")*/
            ("angle,a",
                "Use if the property is an angle and subject to periodic boundary conditions.  ALL PROPERTIES IN THE GIVEN FILE MUST BE ANGLES OR MUST NOT BE ANGLES.  DO NOT MIX ANGLES WITH NON-ANGLES")
            ("bVerbose,v",
                "Print extra information")
            ("bootstrap,b",po::value<double>(&options.nbootstrap)->default_value(0,"Off"),
                "Do bootstrapping --bootstrap times.")
        ;
        try
        {
            std::cout << desc;
            
            po::variables_map vm;
            po::store(po::parse_command_line(argc,argv,desc), vm);
            po::notify(vm);
            
            if (vm.count("help"))
            {
                std::exit(1);
            }
            if (vm.count("bVerbose"))
            {
                options.bVerbose = true;
            }
            if (options.nbootstrap > 0)
            {
                if (options.nbootstrap > 1e15)
                {
                    std::cout << "\n" << std::string(70, '*') << "\n";
                    std::cout << "Cowardly refusing to perform more than 1e15 bootstrap steps." << std::endl;
                    std::cout << std::string(70, '*') << "\n";
                    options.nbootstrap = 1e15;
                }
                options.bootstrap = true;
            }
            if (vm.count("angle"))
            {
                options.isangle = true;
            }
            if (options.datnames.size() != options.datunits.size())
            {
                std::cout << "\nERROR: There are a different number of --dataunit(-u) arguments as --datatype(-d) arguments.\n";
                std::exit(1);
            }
            /*if (vm.count("seed"))
            {
                options.doseed = true;
            }*/
        }
        catch(boost::program_options::required_option& e)
        {
            std::cerr << "\nERROR: " << e.what() << std::endl;
            std::exit(1);
        }
    }
    catch(...)
    {
        std::cerr << "\nException of unknown type(s)\n";
        std::exit(1);
    }
    bw_print_options(options);
    return;
}

void bw_print_options(bw_options &options)
{
    std::cout << std::string(70, '-') << "\n";
    std::cout << "Reading data from " << options.xvglist << std::endl;
    std::cout << "Calculating averages based on data contained in " << options.hdf5file << std::endl;
    std::cout << "Will write data entries as: ";
    for (int i=0; i < options.datnames.size(); i++)
    {
        std::cout << options.datnames[i] << "(" << options.datunits[i] << ") ";
    }
    std::cout << std::endl;
    if (options.isangle)
    {
        std::cout << "Taking into account angle periodicity..." << std::endl;
    }
    if (options.bootstrap)
    {
        std::cout << "Will bootstrap " << options.nbootstrap << " times." << std::endl;
    }
    /*if (options.doseed)
    {
        std::cout << "Using " << options.seed << " as seed" << std::endl;
    }*/
    std::cout << std::string(70, '-') << "\n";
    return;
}
