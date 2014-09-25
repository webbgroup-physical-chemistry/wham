#include "bw_boltzmann_weight.h"

// Read in hdf5 file
// Read in xvg files
// Calculate average

Boltzmann_Weight::Boltzmann_Weight(const bw_options &option) : Interact_H5(option)
{
    options = option;
    bw_read_filelist();
    bw_parse_dat();
    h5_bin_assignments(list);
    bw_calc_prob();
}

void Boltzmann_Weight::calc_average(h5_dat &prob, bw_datfile *list, const int &nexp, bool doHist)
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
    std::vector<double> binprob(nbins,0);
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
    prob.avg = std::vector<double> (nitems,0);
    prob.stdev = std::vector<double> (nitems,0);
    /* Normalize probabilities to 1.0 */
    double sum = 0, invsum = 1;
    for (int i=0; i<nexp; i++)
    {
        for (int j=0; j<(int)list[i].frameN.size(); j+=frame_step)
        {
            frame = list[i].frameN[j];
            if (frame >= prob.span[0] && frame <= prob.span[1] )
            {
                binprob[list[i].bin[frame]] = (double)prob.bin_prob[list[i].bin[frame]] / (double)bincounts[list[i].bin[frame]];
                sum += binprob[list[i].bin[frame]];
                //std::cout << "exp: " << i << ", frame: " << frame << ", bin: " << list[i].bin[frame] << ", counts: " << bincounts[list[i].bin[frame]] << ", prob: " << prob.bin_prob[list[i].bin[frame]] << std::endl;
            }
        }
    }
    invsum = 1/sum;
    std::vector<std::vector<std::vector<double> > > fulldat (nitems);
    // If doing periodic angle average
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
                        x[k] += binprob[list[i].bin[frame]] * invsum * cos(list[i].dat[j][k] * DEG2RAD);
                        y[k] += binprob[list[i].bin[frame]] * invsum * sin(list[i].dat[j][k] * DEG2RAD);
                        if (doHist) {
                            std::vector<double> framedat (2,0);
                            framedat[0] = binprob[list[i].bin[frame]] * invsum;
                            framedat[1] = atan2f(y[k],x[k])*RAD2DEG;
                            fulldat[k].push_back(framedat);
                        }
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
    // If non-periodic average
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
                        total[k] += binprob[list[i].bin[frame]] * invsum * list[i].dat[j][k];
                        if (doHist) {
                            std::vector<double> framedat (2,0);
                            framedat[0] = binprob[list[i].bin[frame]] * invsum;
                            framedat[1] = list[i].dat[j][k];
                            fulldat[k].push_back(framedat);
                        }
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
                        vartotal[k] += binprob[list[i].bin[frame]] * invsum * list[i].dat[j][k] * ( list[i].dat[j][k] - total[k] ) ;
                    }
                }
            }
        }
        for (int i=0;i<nitems;i++)
        {
            prob.avg[i] = total[i];
            prob.stdev[i] = sqrt(vartotal[i]);
            if (doHist) {
                bw_hist(fulldat[i],i);

            }
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
    while (h5_get_dataset(name,probs[n]) == 0)
    {
        if (options.bVerbose)
        {
            std::cout << name << std::endl;
        }
        probs.push_back(placeholder);
        calc_average(probs[n],&list[0],(int)list.size(),false);
        /* set up next set of convergence data */
        n++;
        sprintf(name,"/Ensemble/Conv-%i",n);
    }
    // Then include the dataset over the full range
    h5_get_dataset("/Ensemble",probs[n]);
    if (options.bVerbose)
    {
        std::cout << "/Ensemble" << std::endl;
    }
    calc_average(probs[n],&list[0],(int)list.size(),true);
    if (options.doWrite)
    {
        h5_write_prob(probs);
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
            h5_write_dat(list[i]);
        }
    }
    return;
}

void Boltzmann_Weight::bw_read_dat(bw_datfile &file)
{
    if (options.bVerbose || true)
    {
        std::cout << "Reading experiment " << file.experiment << " from datafile " << file.filename << ".\n";
    }
    xvg_numbered = true;
    std::string line;
    std::vector<double> values;
    double each;
    int i = 0, ncol = 0;
    std::ifstream readfile(file.filename.c_str());
    if (readfile.is_open())
    {
        while (readfile.good())
        {
            getline(readfile,line);
            // Make sure the line is not empty and does not begin with an
            // escape character
            if (not line.empty() && line.substr(0,1) != ";" && line.substr(0,1) != "#" && line.substr(0,1) != "@")
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
                    std::cerr << "\nERROR: The number of columns (" << values.size() << ") in " << file.filename << " are not the same number of items give with --datatype (" << options.datnames.size() << ").  The offending line is:\n" << std::endl;
                    std::cerr << line << std::endl << std::endl;;
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

void Boltzmann_Weight::bw_hist(const std::vector<std::vector<double> > &dat, int &index)
{
    int frames = dat.size();
    double min = dat[0][1], max = dat[0][1];
    for (int i = 0; i<frames; i++) {
        if (dat[i][1] > max) { max = dat[i][1]; }
        if (dat[i][1] < min) { min = dat[i][1]; }
    }
    max += MAX(abs(max*0.01),0.1);
    min -= MAX(abs(min*0.01),0.1);
    double range = max - min;
    double stepsize = range / options.histbins;
    std::vector<std::vector<double> > histogram (options.histbins, std::vector<double> (2,0.));
    for (int i = 0; i<frames; i++) {
        int n = std::floor((dat[i][1] - min) / stepsize);
        if (n < 0) {
            
            std::cerr << "Binning index too LOW: "<< i << " " << (dat[i][1] - min) / stepsize << " " << n << " " << dat[i][1] << " " << min << " " << max << std::endl;
            std::exit(1);
        }
        else if (n >= options.histbins) {
            std::cerr << "Binning index too HIGH: "<< i << " " << (dat[i][1] - min) / stepsize << " " << n << " " << dat[i][1] << " " << min << " " << max << std::endl;
            std::exit(1);
        }
        double leftedge = n*stepsize+min;
        double rightedge = (n+1)*stepsize+min;
        if (dat[i][1] < leftedge ) {
            std::cerr << "\nError in binning (smaller than this bin edge)\n";
            std::cerr << n << " " << leftedge << " " << dat[i][1] << " " << rightedge << std::endl;
            std::exit(1);
        }
        else if (dat[i][1] >= rightedge ) {
            std::cerr << "\nError in binning (larger than next bin edge)\n";
            std::cerr << n << " " << leftedge << " " << dat[i][1] << " " << rightedge << std::endl;
            std::exit(1);
        }
        histogram[n][1] += dat[i][0];
    }
    // Assign all left edges; did not do in previous loop b/c some left edges
    // may not be visited.
    for (int i=0;i<options.histbins;i++){
        histogram[i][0] = i * stepsize + min;
    }
    h5_write_hist(histogram,index);
    return;
}
 