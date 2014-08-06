#include "parse_options.h"


void wham_options(int argc, char *argv[], t_options &options)
{
    namespace po = boost::program_options;
    const char *description =
    {
        "\tThis is a multidimensional Weighted-Histrogram Analysis program \nwritting in C++.  For more information about the algorithm employed, \nrefer to 'BEFEUS: Bayesian Estimation of Free Energy for Umbrella \nSampling' by Dr. Daniel L. Ensign. \n\tCurrently, this build can only handle up to 3 degrees of freedom, \nbut has only been tested up to 2. \n\tNOTE 1: Be sure to make sure the file list format correctly adheres \nto what is listed under the --filelist option information.  The trajectory # \nneeds to be unique for each individual simulation.  This parameter is \nused to group different degrees of freedom to the same biased \ntrajectory.  The numbering can start at 0 -or- 1. \n\tNOTE 2:  When using the --first_frame and --last_frame flags, the \nselection will show up at the begining and the end of the run, but the number \nof frames in each trajectory will still be listed as the total number of \nframes for the trajectory, rather than however many were used. \n\tNOTE 3: When assigning bin widths and coordinate ranges, use multiple \ninvocations of -w or include the values as a string within quotations (the \norder is the same order as the dimensions are listed in the --filelist file). \n\tNOTE 4: Make sure the unit (X) on --min/max_coordinate and --bin_width \nmatch the coordinate bias value, flat range, and force constant in the \n--filelist input.\n\tNOTE 4: Assuming every trajectory has the same number of frames.  I can fix this assumption, but for now it's easier not to. \n\nOptions"
    };
    std::vector<std::string> binws;
    std::vector<std::string> xmin;
    std::vector<std::string> xmax;
    options.bVerbose = false;
    options.doConv = false;
    options.nexceed = 0;
    try
    {
        po::options_description desc(description);
        desc.add_options()
            ("help,h","")
            ("filelist,f",po::value<std::string>(&options.filelist)->required(),
                "File containting a list of trajectories and information about the umbrella potentials applied.\nFromat: <trajectory #> <coordinate file> <coordinate bias value (X)> <flat range (X)> <force constant (kJ mol^-1 X^-2)> <temperature (K)>")
            ("out,o",po::value<std::string>(&options.outname),
                "Output file name")
            ("first_frame,b",po::value<int>(&options.f0)->default_value(0),
                "First frame to read")
            ("last_frame,e",po::value<int>(&options.fN)->default_value(-1,"last frame"),
                "Last frame to read, where <-e> is an integer")
            ("temperature,t",po::value<float>(&options.t)->default_value(300),
                "Temperature (K)")
            ("min_coordinate,n",po::value<std::vector<std::string> >(&xmin),
                "Minimum value of reaction coordinate. Default: -180 (X)")
            ("max_coordinate,x",po::value<std::vector<std::string> >(&xmax),
                "Maximum value of reaction coordinate. Default:  180 (X)")
            ("bin_width,w",po::value<std::vector<std::string> >(&binws),
                "Bin width. Default: 5 (X)")
            ("dof,d",po::value<int>(&options.ndof)->default_value(1),
                "Number of degrees of freedom being biased")
            ("tolerance,l",po::value<float>(&options.tol)->default_value(1e-12,"1e-12"),
                "Convergence tolerance for monte carlo")
            ("iterations,i",po::value<float>(&options.iter)->default_value(5e5,"5e5"),
                "Maximum number of iterations")
            ("convergence,c",po::value<int>(&options.convStep)->default_value(0,"Off"),
                "Perform WHAM every <-c> steps, where <-c> is an integer")
            ("bVerbose,v",
                "Print extra information")
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
            if (options.convStep > 0)
            {
                options.doConv = true;
            }
            /* The default output name is the input file name suffixed with the appropriate extensions */
            if (!vm.count("out"))
            {
                options.outname = options.filelist;
            }
            /* Make sure the outfile suffix is correct */
            if ( options.outname.find(".h5") == std::string::npos)
            {
                //if (options.outname.back() == '.' ) // -std=c++11, which isn't working with boost and/or hdf5
                if (*options.outname.rbegin() == '.' )
                    {
                        options.outname.append("h5");
                    }
                else
                {
                    options.outname.append(".h5");
                }
            }
            /* Allow for different bin widths and coordinate ranges in each dimension */
            multidimensional_option(binws,options,5.,"--bin_width",options.b);
            multidimensional_option(xmin,options,-180.,"--min_coordinate",options.x0);
            multidimensional_option(xmax,options,180.,"--max_coordinate",options.xN);
             
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

    print_options(options);
    return;
}

void multidimensional_option(std::vector<std::string> option,
                             t_options &options,
                             float defaultX,
                             std::string valuename,
                             std::vector<float> &X)
{
    if (option.size() > 0)
    {
        for (size_t i=0; i<option.size(); i++)
        {
            try
            {
                float bw = boost::lexical_cast<float>(option[i]);
                X.push_back(bw);
            }
            catch (boost::bad_lexical_cast const&)
            {
                std::cerr << "You entered '" << option[i] << "' as a bin width.  This should be a number" << std::endl;
                std::exit(1);
            }
        }
    }
    else
    {
        X.push_back(defaultX);
    }
    if ((int)X.size() < options.ndof)
    {
        std::cout << "\nThe number of degrees of freedom (" << options.ndof <<") and the number of " << valuename << "s (" << X.size() << ") do not match.  Using the first " << valuename << ", " << X[0] << ", for each dimension.\n";
        for (int i=1; i<options.ndof; i++)
        {
            X.push_back(X[0]);
        }
    }
    else if ((int)X.size() > options.ndof)
    {
        std::cerr << "\nToo many " << valuename << "s given.  std::exiting to avoid assumptions." << std::endl;
        std::exit(1);
    }
    
}

void print_options(t_options &options)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "Will read <" << options.filelist << "> from frame " << options.f0 << " to " << options.fN << "." << std::endl;
    std::cout << "Temperature = " << options.t << " K" << std::endl;
    std::cout << "Thermal Energy = " << KB*options.t << " kJ/mol" << std::endl;
    std::cout << "Biasing coordinate range(s) : Bin width(s) = \n";
    for (int i=0; i<options.ndof; i++)
    {
        std::cout << "\t\t" << options.x0[i] << " to " << options.xN[i] << " : " << options.b[i] << std::endl;;
    }
    std::cout << "Degrees of freedom = " << options.ndof << std::endl;
    std::cout << "Will iterate up to " << options.iter << " times or until converged to " << options.tol <<"."<< std::endl;
    std::cout << "Will write to " << options.outname << std::endl;
#ifdef _OPENMP
    std::cout << "NTHREADS = " << omp_get_max_threads() << std::endl;
#endif
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}
