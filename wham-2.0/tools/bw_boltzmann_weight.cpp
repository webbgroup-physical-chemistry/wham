#include "bw_boltzmann_weight.h"

// Read in hdf5 file
// Read in xvg files
// Calculate average

void Boltzmann_Weight::bw_init(bw_options option)
{
    options = option;
    h5file.h5_init(options);
    bw_read_filelist();
    bw_parse_dat();
    h5file.h5_bin_assignments(list);
    bw_calc_prob();
    return;
}

void Boltzmann_Weight::calc_average(h5_dat &prob, bw_datfile *list, int nexp)
{
    /* 
     * To easily play around with and determine if skipping frames uniformly returns
     * the same answer.  It does not.
     */
    int frame_step = options.frameStep;
    if (options.bVerbose)
    {
        std::cout << "Frames " << prob.span[0]+1 << "-" << prob.span[1] << ".\n";
    }
    /* count how many times each bin is visited in the dataset */
    int nbins = prob.bin_prob.size(),frame;
    std::vector<int> bincounts(nbins,0);
    for (int i=0; i<nexp; i++)
    {
        for (int j=0; j<(int)list[i].frameN.size(); j+=frame_step)
        {
            frame = list[i].frameN[j];
            if (frame >= prob.span[0] && frame <= prob.span[1])
            {
                bincounts[list[i].bin[frame]]++;
            }
        }
    }
    int nitems = list[0].dat[0].size();
    std::vector<float> ph(nitems,0);
    prob.avg = ph;
    prob.stdev = ph;
    /* Normalize probabilities to 1.0 */
    double sum = 0, invsum = 1;
    for (int i=0; i<nexp;i++)
    {
        for (int j=0; j<(int)list[i].frameN.size(); j+=frame_step)
        {
            frame = list[i].frameN[j];
            if (frame >= prob.span[0] && frame <= prob.span[1] )
            {
                sum += (double)prob.bin_prob[list[i].bin[frame]] / (double)bincounts[list[i].bin[frame]];
                //std::cout << "dat: " << prob.bin_prob[list[i].bin[frame]] << " " << frame << " " << bincounts[list[i].bin[frame]] << " " << prob.bin_prob[list[i].bin[frame]] / bincounts[list[i].bin[frame]] << std::endl;
                //std::cout << "dat: " << i << ", frame: " << frame << ", bin: " << list[i].bin[frame] << ", bin prob: " << prob.bin_prob[list[i].bin[frame]] << ", bincounts: " << bincounts[list[i].bin[frame]] << ", frame prob: " << prob.bin_prob[list[i].bin[frame]] / bincounts[list[i].bin[frame]] << " " << list[i].dat[j][k] << std::endl;
            }
        }
    }
    invsum = 1/sum;
    
    if (options.isangle)
    {
        std::vector<double> x(nitems,0), y(nitems,0);
        for (int i=0; i<nexp;i++)
        {
            for (int j=0; j<(int)list[i].frameN.size(); j+=frame_step)
            {
                frame = list[i].frameN[j];
                if (frame >= prob.span[0] && frame <= prob.span[1])
                {
                    for (int k=0; k<nitems; k++)
                    {
                        x[k] += prob.bin_prob[list[i].bin[frame]] / bincounts[list[i].bin[frame]] * invsum * cos(list[i].dat[j][k] * DEG2RAD);
                        y[k] += prob.bin_prob[list[i].bin[frame]] / bincounts[list[i].bin[frame]] * invsum * sin(list[i].dat[j][k] * DEG2RAD);
                    }
                }
            }
        }
        for (int i=0; i<nitems; i++)
        {
            prob.avg[i] = atan2f(y[i],x[i])*RAD2DEG;
            prob.stdev[i] = sqrt( (1-sqrt(x[i]*x[i] + y[i]*y[i])) * 360);
        }
    }
    else
    {
        std::vector<double> total(nitems,0),vartotal(nitems,0);
        for (int i=0; i<nexp;i++)
        {
            for (int j=0; j<(int)list[i].frameN.size(); j+=frame_step)
            {
                frame = list[i].frameN[j];
                if (frame >= prob.span[0] && frame <= prob.span[1])
                {
                    for (int k=0; k<nitems; k++)
                    {
                        total[k] += prob.bin_prob[list[i].bin[frame]] / bincounts[list[i].bin[frame]] * invsum * list[i].dat[j][k];
                        //std::cout << "dat: " << i << " " << frame << " " << list[i].bin[frame] << " " << prob.bin_prob[list[i].bin[frame]] << " " << bincounts[list[i].bin[frame]] << " " << sum << " " << invsum << " " << list[i].dat[j][k] << std::endl;
                    }
                }
            }
        }
        for (int i=0; i<nexp;i++)
        {
            for (int j=0; j<(int)list[i].frameN.size(); j+=frame_step)
            {
                frame = list[i].frameN[j];
                if (frame >= prob.span[0] && frame <= prob.span[1])
                {
                    for (int k=0; k<nitems; k++)
                    {
                        vartotal[k] += prob.bin_prob[list[i].bin[frame]] / bincounts[list[i].bin[frame]] * invsum * list[i].dat[j][k] * ( list[i].dat[j][k] - total[k] ) ;
                    }
                }
            }
        }
        for (int i=0;i<nitems;i++)
        {
            prob.avg[i] = total[i];
            prob.stdev[i] = sqrt(vartotal[i]);
        }
    }
    return;
}

void Boltzmann_Weight::bw_calc_prob()
{
    // Loop through all convergence data sets until we no longer find one
    std::vector<h5_dat> probs;
    h5_dat placeholder;
    char name[1024];
    int n = 0;
    sprintf(name,"/Ensemble/Conv-%i",n);
    probs.push_back(placeholder);
    while (h5file.h5_get_dataset(name,probs[n]) == 0)
    {
        if (options.bVerbose)
        {
            std::cout << name << std::endl;
        }
        probs.push_back(placeholder);
        calc_average(probs[n],&list[0],(int)list.size());
        /* set up next set of convergence data */
        n++;
        sprintf(name,"/Ensemble/Conv-%i",n);
    }
    // Then include the dataset over the full range
    h5file.h5_get_dataset("/Ensemble",probs[n]);
    if (options.bVerbose)
    {
        std::cout << "/Ensemble" << std::endl;
    }
    calc_average(probs[n],&list[0],(int)list.size());
    if (options.doWrite)
    {
        h5file.h5_write_prob(probs);
    }
    // Extra stuff for some scripts
    std::cout << "Final Averages:\n";
    for (int i=0; i<(int)probs[n].avg.size(); i++)
        std::cout << "Average: " << probs[n].avg[i] << ", STDEV: " << probs[n].stdev[i] << ", DATASET: " << options.datnames[i] << ", UNITS: " << options.datunits[i] << std::endl;
    for (int i=0; i<(int)probs[n].avg.size(); i++)
    {
        if (i != (int)probs[n].avg.size() -1)
        {
            std::cout << probs[n].avg[i] << " ";
        }
        else
        {
            std::cout << probs[n].avg[i] << std::endl;
        }
    }
}

void Boltzmann_Weight::bw_parse_dat()
{
    for (int i=0; i<(int)list.size(); i++)
    {
        bw_read_dat(list[i]);
        if (options.bVerbose)
        {
            std::cout << "\t" << list[i].dat.size() << " frames  with " << list[i].dat[0].size() << " columns in each row. " << std::endl;
        }
        if (options.doWrite)
        {
            h5file.h5_write_dat(list[i]);
        }
    }
    return;
}

void Boltzmann_Weight::bw_read_dat(bw_datfile &file)
{
    if (options.bVerbose)
    {
        std::cout << "Reading experiment " << file.experiment << " from datafile " << file.filename << ".\n";
    }
    xvg_numbered = true;
    std::string line;
    std::vector<float> values;
    float each;
    int i = 0, ncol = 0;
    std::ifstream readfile(file.filename.c_str());
    if (readfile.is_open())
    {
        while (readfile.good())
        {
            getline(readfile,line);
            if (not line.empty() && line.substr(0,1) != ";" && line.substr(0,1) != "#")
            {
                std::stringstream linestream(line);
                while (linestream >> each)
                {
                    values.push_back(each);
                }
                if (values.size() == options.datnames.size()+1 && xvg_numbered)
                {
                    // Frame number is included
                    if ( values[0] >= 0 )
                    {
                        file.frameN.push_back(values[0]);
                        values.erase(values.begin() + 0);
                    }
                    else
                    {
                        std::cerr << "\nERROR: The frame number in " << file.filename << ", line " << i << ", is negative.\n" << line << std::endl;
                        std::exit(1);
                    }
                }
                else if (values.size() == options.datnames.size())
                {
                    // Frame number is excluded
                    if (xvg_numbered)
                    {
                        std::cout << "\nWARNING: Frame numbers not specified in " << file.filename << ".  Assuming all frames are included, starting from frame zero." << std::endl;
                        xvg_numbered = false;
                    }
                    file.frameN.push_back(i);
                }
                else
                {
                    // Who knows...
                    std::cerr << "\nERROR: The number of columns (" << values.size() << ") in " << file.filename << " are not the same number of items give with --datatype (" << options.datnames.size() << ")." << std::endl;
                    if ((int)values.size() > 1)
                    {
                        std::cerr << "To include multiple entries to --datatype, invoke the option multiple times in the same order as the columns from left-to-right." << std::endl;
                    }
                    std::exit(1);
                }
                if (i == 0)
                {
                    ncol = values.size();
                }
                else
                {
                    if ( (int)values.size() != ncol )
                    {
                        std::cerr << "\nERROR: The number of columns (" << values.size() << ") in " << file.filename << " are not the same in row " << i << " as the first row." << std::endl;;
                        std::exit(1);
                    }
                }
                file.dat.push_back(values);
                // Everything is good here...
                //std::cout << i << " " << file.frameN[i] << " " << file.dat[i][0] << " (" << file.frameN.size() << ")" << std::endl;
                values.clear();
                values.shrink_to_fit();
                i++;
            }

        }
    }
    else
    {
        std::cerr << "\nError: Cannot open " << file.filename << ". Exiting..." << std::endl;
        std::exit(1);
    }
    if ( i < 1 )
    {
        std::cerr << "\nError: " << file.filename << " is empty.  Reconsider the contents of " << options.xvglist << ".\n";
        std::exit(1);
    }
}

void Boltzmann_Weight::bw_read_filelist()
{
    std::string line;
    std::vector<std::string> values;
    bw_datfile item;
    std::string name;
    std::ifstream file(options.xvglist.c_str());
    int n = 0;
    if (file.is_open())
    {
        if (options.bVerbose)
        {
            std::cout << "Reading " << options.xvglist << ".\n";
        }
        while (file.good())
        {
            getline( file, line);
            if (not line.empty())
            {
                std::stringstream linestream(line);
                while (linestream >> name)
                {
                    values.push_back(name);
                    if ((int)values.size() > 2)
                    {
                        std::cerr << "\nERROR:  Too many items on line " << n << " in " << options.xvglist << ".  Should be only:\n\tExperiment# datafile\n" << std::endl;
                        std::exit(1);
                    }
                }
                int number;
                try
                {
                    number = boost::lexical_cast<int>(values[0]);
                    item.experiment = number;
                    item.filename = values[1];
                }
                catch(boost::bad_lexical_cast& e)
                {
                    std::cerr << "\nERROR:  Expected an integer argument on line " << n << " indicating the experiment number.  Instead, recieved " << values[0] << std::endl;
                    std::exit(1);
                }
                list.push_back(item);
                values.clear();
                values.shrink_to_fit();
                n++;
            }
        }
    }
    else
    {
        std::cerr << "Error: Cannot open "<< options.xvglist << ". Exiting..." << std::endl;
        std::exit(1);
    }
    return;
}
