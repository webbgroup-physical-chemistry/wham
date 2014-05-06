#include "wham_befeus.h"

void TorsionExperiment::read_filelist_and_setup( std::string line, t_options parser_options, int previous_experiment, int previous_dof)
{
    /* Assign the parser options */
    options = parser_options;
    
    std::string item;
    std::stringstream linestream(line);
    int nitems = 0;
    while (linestream >> item)
    {
        nitems++;
    }
    if (nitems != 6)
    {
        std::cerr << "\nError! Input file incorrectly formatted.  Correct format is:\n<experiment number> <xvg path> <biased position> <zero potential +/- range> <force constant> <temperature>" << std::endl;
        std::exit(1);
    }
    /* Assign values form the file list to the experiment */
    std::stringstream relinestream(line);
    relinestream >> experiment >> xvgfile >> phi0 >> deltaphi >> fc >> t;
    beta = 1/(KB*t);
    fc *= DEG2RAD * DEG2RAD; //convert from kJ/mol/radian^2 to kJ/mol/degree^2
    /* If the experiment number is the same as the previous experiment number, it's an additional degree of freedom.  Otherwise, we reset back to the zeroth degree of freedom */
    if (experiment == previous_experiment)
    {
        dof = previous_dof + 1;
    }
    else
    {
        dof = 0;
    }
    if (dof >= options.ndof)
    {
        std::cerr << "\nERROR: Found " << dof+1 << " degrees of freedom when expecting " << options.ndof << " degrees of fredom.\n" << std::endl;
        std::exit(1);
    }
    /* Now we can assign the ranges and bin sizes to match the specific DoF */
    begin = options.x0[dof];
    end = options.xN[dof];
    bin_size = options.b[dof];
    
    if (options.bVerbose)
    {
        std::cout << "\nReading data file: " << xvgfile << " from trajectory " << experiment << ", degree of freedom #" << dof+1 << std::endl;
        std::cout << "Biasing coordinate centered on " << phi0 << " +/- " << deltaphi << std::endl;
        std::cout << "Data ranges from " << begin << " to " << end << ".  Will look at bins " << bin_size << " units wide" << std::endl;
        std::cout << "Force constant = " << fc << " kJ/mol/degree^2 at temperature " << beta << " mol*kJ^-1" << std::endl;
    }
    precision_error = false;
    set_nbins();
    read_trajectory();
    build_histogram();
    assign_omega();
    weight_dihedral_restraint_potentials();
    if (precision_error)
    {
        std::cerr << "\nWARNING!  Numeric precision, " << std::numeric_limits<double>::min() << ", has been exceeded.  " << std::endl;
    }
    return;
}

void TorsionExperiment::read_trajectory()
{
    std::string line;
    std::ifstream file(xvgfile.c_str());
    if (file.is_open())
    {
        while (file.good())
        {
            getline( file, line);
            if (not line.empty())
            {
                /* Since this was writen with gromacs in mind, the input xvgfile has 2 columns -- time coordinate -- and here we are checking that format.  It may also be convenient to to have a single-column data file, so we check these 2 cases.  Otherwise, it will fail. */
                std::stringstream check_line(line);
                std::string item;
                int nitems = 0;
                float time = 0, value = 0;
                while (check_line >> item)
                {
                    nitems++;
                }
                std::stringstream linestream(line);
                switch (nitems)
                {
                    case 1 :
                        linestream >> value;
                        coordinate.push_back(value);
                        break;
                    case 2 :
                        linestream >> time >> value;
                        coordinate.push_back(value);
                        break;
                    default :
                        std::cerr << "\nError!  There are " << nitems << " columns in " << xvgfile << ".  Was expecting 1 (data) or 2 (time data) columns" << std::endl;
                        std::exit(1);
                }
            }
        }
    }
    else if (!file)
    {
        std::cerr << "\nError reading " << xvgfile << std::endl;
        std::exit(1);
    }
    file.close();
    return;
}

void TorsionExperiment::set_nbins()
{
    float range = end - begin;
    float resid = remainder(range,bin_size);
    if (resid < 1e-9)
    {
        nbins = boost::numeric_cast<int>(range/bin_size);
    }
    else
    {
        std::cerr << "\nWARNING: There is not an integer number of bins using the user assigned bin width (" << bin_size << ") over the range " << begin << " - " << end << ".  Using a default of 50 bins with a width of " << range/50 << " units." << std::endl;
        bin_size = range / 50;
        nbins = 50;
    }
    return;
}
void TorsionExperiment::build_histogram()
{
    int ncoords = coordinate.size();
    std::vector<t_bin> fb(ncoords);
    frame_bins = fb;
    int bin = 0;
    for (int i=begin; i<end; i+=bin_size)
    {
        bincounts.push_back(0);
        for (int j=0; j<ncoords; j++)
        {
            if (coordinate[j] >= begin && coordinate[j] <= end)
            {
                if ( i <= coordinate[j] && coordinate[j] < i+bin_size)
                {
                    bincounts[bin]++;
                    frame_bins[j].bin1D = bin;
                }
            }
            else
            {
                std::cerr << "\nERROR: frame " << j << " has a value of " << coordinate[j] << " which is outside the range of " << begin << " - " << end << std::endl;
                std::exit(1);
            }
        }
        bin++;
    }
    int nbins = bincounts.size();
    if (options.bVerbose)
    {
        std::cout << "1-Dimensional bin counts: [ ";
        for (int i=0; i<nbins; i++)
        {
            std::cout << bincounts[i] << " ";
        }
        std::cout << "]" << std::endl;
    }
    return;
}

void TorsionExperiment::weight_dihedral_restraint_potentials()
{
    int ncoords = coordinate.size();
    for (int i=0; i<ncoords; i++)
    {
        potential.push_back(weight_dihedral_restraint_potential(coordinate[i]));
    }
    return;
}

float TorsionExperiment::weight_dihedral_restraint_potential(float phi)
{
    float diffphi = delta_angle(phi0,phi);
    float ddp;
    if ( diffphi > deltaphi )
    {
        ddp = delta_angle( diffphi, deltaphi);
    }
    else
    {
        ddp = 0;
    }
    return exp(-beta * 0.5 * fc * ddp * ddp);
}

double TorsionExperiment::simpson_integration(float a, float b, int steps)
{
    double integrand = 0;
    double width = ((double)b - (double)a)/(double)steps;
    double fa, fb, fab;
    for (double i=(double)a; i<=(double)b; i+=width)
    {
        fa = weight_dihedral_restraint_potential(i);
        fb = weight_dihedral_restraint_potential(i + width);
        fab = weight_dihedral_restraint_potential((i+i+width)/2);
        integrand += width/6 * (fa + 4 * fab + fb);
    }
    return integrand;
}

double TorsionExperiment::trapezoidal_integration( float a, float b, int steps)
{
    double integrand = 0;
    double width = ((double)b - (double)a)/(double)steps;
    double fa, fb;
    for (float i=(double)a; i<=(double)b; i+=width)
    {
        fa = weight_dihedral_restraint_potential(i);
        fb = weight_dihedral_restraint_potential(i+width);
        integrand += 0.5 * width * (fa+fb);
    }
    return integrand;
}

void TorsionExperiment::assign_omega()
{
    for (int i=0; i<nbins; i++)
    {
        float a = i * bin_size + begin;
        float b = (i+1) * bin_size + begin;
        omega.push_back(simpson_integration(a,b,1e2));
        if ( omega[i] <= std::numeric_limits<double>::min() )
        {
            precision_error = true;
        }
    }
    return;
}

void TorsionExperiment::assign_NDbin(int frame, int bin)
{
    frame_bins[frame].binND = bin;
    return;
}

int TorsionExperiment::Experiment()
{
    return experiment;
}

int TorsionExperiment::DoF()
{
    return dof;
}

int TorsionExperiment::nBins()
{
    return nbins;
}

float TorsionExperiment::T()
{
    return t;
}

float TorsionExperiment::Phi0()
{
    return phi0;
}

float TorsionExperiment::DPhi()
{
    return deltaphi;
}

float TorsionExperiment::FC()
{
    return fc;
}

std::vector<float> TorsionExperiment::Potential()
{
    return potential;
}

std::vector<double> TorsionExperiment::Omega()
{
    return omega;
}

std::vector<float> TorsionExperiment::Coordinate()
{
    return coordinate;
}

std::vector<int> TorsionExperiment::BinCounts()
{
    return bincounts;
}

std::vector<t_bin> TorsionExperiment::FrameBins()
{
    return frame_bins;
}