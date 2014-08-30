#include "wham.h"

/*
void cpp_wham_conv(WHAM wham, t_options options)
{
    std::cerr << "\nDoes not work correctly.  Exiting for your sanity...\n";
    std::exit(1);
    const int nConv = ( wham.get_num_frames() - options.f0 ) /options.convStep - 1;
    std::vector<WHAM> converge(nConv);
    
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int n=0 ; n<nConv; n++)
    {
        converge[n] = wham;
        int f0 = converge[n].get_options_f0();
        int fN = f0 + n * options.convStep;
        converge[n].set_options_fN(fN);
        std::clock_t start;
        if (converge[n].get_options_f0() != converge[n].get_options_fN() && converge[n].get_options_fN() != converge[n].get_final_frame())
        {
            start = std::clock();
            converge[n].cpp_wham_initialize_vectors();
            converge[n].cpp_wham_count_bins();
            converge[n].cpp_wham_mapOmega();
            converge[n].cpp_wham_doWHAM();
            //Write_HDF5 outhdf5 = wham.get_outhdf5();
            //char groupname[1024];
            //sprintf(groupname,"/Ensemble/Conv-%d",n);
            //std::vector<float> prob = converge[n].get_opt_prob();
            //std::vector<float> pmf  = converge[n].get_opt_pmf();
            //t_wham wham_args = converge[n].get_wham_args();
            //std::vector<t_map> map = converge[n].get_map();
            //outhdf5.write_wham(prob,pmf,wham_args.counts,map,groupname,f0,fN);
            if (options.bVerbose || true)
            {
#ifdef _OPENMP
                fprintf(stdout,"THREAD-%d: Completed frames %d to %d in %.1f ms\n",omp_get_thread_num(),f0,fN,(std::clock()-start)/(float)(CLOCKS_PER_SEC/1000));
#else
                fprintf(stdout,"Completed frames %d to %d in %.1f ms\n",f0,fN,(std::clock()-start)/(float)(CLOCKS_PER_SEC/1000));
#endif
            }
        }
    }
    Write_HDF5 outhdf5 = wham.get_outhdf5();
    for (int n=0; n<nConv; n++)
    {
        char groupname[1024];
        sprintf(groupname,"/Ensemble/Conv-%d",n);
        std::vector<float> prob = converge[n].get_opt_prob();
        std::vector<float> pmf  = converge[n].get_opt_pmf();
        t_wham wham_args = converge[n].get_wham_args();
        std::vector<t_map> map = converge[n].get_map();
        int f0 = converge[n].get_options_f0();
        int fN = converge[n].get_options_fN();
        outhdf5.write_wham(prob,pmf,wham_args.counts,map,groupname,f0,fN);
    }
    return;
}
*/

void WHAM::cpp_wham_conv()
{
    if (options.doConv)
    {
        char groupname[1024];
        int n = 0;
        int final_frame = last_frame;
        std::clock_t start;
        for (int i=options.f0; i<nframes; i += options.convStep)
        {
            options.fN = i;
            if (options.f0 != options.fN && options.fN != final_frame)
            {
                start = std::clock();
                std::cout << "\nWHAM over frames " << options.f0 << " to " << i;
                cpp_wham_initialize_vectors();
                cpp_wham_count_bins();
                cpp_wham_mapOmega();
                cpp_wham_doWHAM();
                sprintf(groupname,"/Ensemble/Conv-%i",n);
                outhdf5.write_wham(opt_prob, opt_pmf, wham_args.counts, map, groupname, options.f0, options.fN);
                if (options.bVerbose || true)
                {
                    std::cout << "Completed frames " << options.f0 << "-" << i << " in " << (std::clock()-start)/(float)(CLOCKS_PER_SEC/1000) << " ms\n";
                }
                n++;
            }
        }
    }
}

void WHAM::cpp_wham_init(t_options option)
{
    options = option;
    cpp_wham_read_experiments();
    cpp_wham_find_nstates();
    cpp_wham_initialize_vectors();
    cpp_wham_group_ndim();
    cpp_wham_make_ND_map();
    cpp_wham_assign_bins();
    cpp_wham_count_bins();
    cpp_wham_mapOmega();
    cpp_wham_doWHAM();
    outhdf5.init(options, wham_args, nexp, nstates, opt_prob, opt_pmf);
    outhdf5.write_trj(group_traj);
    outhdf5.write_wham(opt_prob, opt_pmf, wham_args.counts, map, "/Ensemble", options.f0, last_frame);
}

void WHAM::cpp_wham_read_experiments()
{
    /* Initialize vector of experiments */
    experiment = std::vector<TorsionExperiment> (options.ntraj);
    std::string line;
    std::ifstream file(options.filelist.c_str());
    int i = 0;
    first_experiment = 0;
    int previous_experiment = -1;
    int previous_dof = -1;
    if (file.is_open())
    {
        if (options.bVerbose)
        {
            std::cout << "Opening " << options.filelist << std::endl;
        }
        while (file.good())
        {
            getline( file, line );
            if (not line.empty())
            {
                experiment[i].read_filelist( line, options, previous_experiment, previous_dof );
                previous_experiment = experiment[i].Experiment();
                previous_dof = experiment[i].DoF();
#ifndef _OPENMP
                experiment[i].setup(options);
#endif
                i++;
            }
        }
    }
    else if (!file)
    {
        std::cerr << "\nError opening " << options.filelist << "." << std::endl;
        std::exit(1);
    }
    file.close();

    first_experiment = experiment[0].Experiment();

#ifdef _OPENMP
    #pragma omp parallel for num_threads(omp_get_max_threads()) 
    for (int n=0;n<i;n++)
    {
        experiment[n].setup(options);
    }
#endif

    return;
}

void WHAM::cpp_wham_find_nstates()
{
    nstates = 1;
    for (int i=0; i<options.ndof; i++)
    {
        for (int j=0; j<options.ntraj; j++)
        {
            if ( experiment[j].DoF() == i )
            {
                nstates *= experiment[j].nBins();
                break;
            }
        }
    }
    invnstates = 1./nstates;
    nexp = options.ntraj/options.ndof;
    nframes = 0;
    for (int i=0; i<nexp; i++)
    {
        int trajnframes = experiment[i].Coordinate().size();
        if (nframes == 0 || nframes > trajnframes)
        {
            nframes = trajnframes;
        }
    }
}

void WHAM::cpp_wham_initialize_vectors()
{
    // Clearing all vectors, should be useful for convergence
    wham_args.omegas.clear();
    wham_args.sample.clear();
    wham_args.counts.clear();
    wham_args.t.clear();
    
    wham_args.sample = std::vector<int> (nexp,0);
    wham_args.counts = std::vector<int> (nstates,0);
    wham_args.t = std::vector<float> (nexp,300);
    wham_args.omegas = std::vector<std::vector<double> > (nexp,std::vector<double> (nexp,1));

    return;
}

void WHAM::cpp_wham_group_ndim()
{
    std::string grammar = "";
    if (options.ndof > 1)
    {
        grammar = "s";
    }
    
    bool first_traj = false;
    std::vector<t_experiment> group_traj2(nexp);
    // pieces, [dofM]
    
    int nn = 0;
    for (int i=first_experiment; i<(first_experiment+nexp); i++)
    {
        int ndof = 0;
        int expn,expdof;

        for (int j=0; j<options.ntraj; j++)
        {
            expn = experiment[j].Experiment();
            if (expn == i )
            {
                first_traj = false;
                // Store these since we use them a couple times
                expdof = experiment[j].DoF();
                std::vector<float> expcoord = experiment[j].Coordinate();
                std::vector<float> exppot = experiment[j].Potential();
                std::vector<t_bin> expbin = experiment[j].FrameBins();
                std::vector<double> expomega = experiment[j].Omega();
                std::vector<int> expbincount = experiment[j].BinCounts();
                int ncoords = expcoord.size();
                int nomegas = expomega.size();
                int expnbins = experiment[j].nBins();

                if (expdof == 0)
                {
                    first_traj = true;
                }
                // Initialize the vectors for the first trajectory in a DoF only
                if (first_traj)
                {
                    group_traj2[nn].nbins = std::vector<int> (options.ndof);
                    group_traj2[nn].phi0 = std::vector<float> (options.ndof);
                    group_traj2[nn].deltaphi = std::vector<float> (options.ndof);
                    group_traj2[nn].fc = std::vector<float> (options.ndof);

                    group_traj2[nn].experiment = i;
                }
                group_traj2[nn].nbins[expdof] = expnbins;
                group_traj2[nn].phi0[expdof] = experiment[j].Phi0();
                group_traj2[nn].deltaphi[expdof] = experiment[j].DPhi();
                group_traj2[nn].fc[expdof] = experiment[j].FC();
                if (group_traj2[nn].t > 0 && experiment[j].T() != group_traj2[nn].t ) {
                    std::cerr << "\nERROR! Temperature don't match for experiment " << i << "." << std::endl;
                }
                group_traj2[nn].t = experiment[j].T();
                if (first_traj)
                {
                    group_traj2[nn].potential = std::vector<std::vector<float> >(ncoords, std::vector<float> (options.ndof,0));
                    group_traj2[nn].coordinate = std::vector<std::vector<float> >(ncoords, std::vector<float> (options.ndof,0));
                    group_traj2[nn].frame_bin = std::vector<std::vector<t_bin> > (ncoords, std::vector<t_bin> (options.ndof));
                    
                    group_traj2[nn].omega = std::vector<std::vector<double> >(expnbins, std::vector<double> (options.ndof,1));
                    group_traj2[nn].bincounts = std::vector<std::vector<int> >(expnbins, std::vector<int> (options.ndof,0));
                }
                for (int k=0; k<ncoords; k++)
                {
                    group_traj2[nn].potential[k][expdof] = exppot[k];
                    group_traj2[nn].coordinate[k][expdof] = expcoord[k];
                    group_traj2[nn].frame_bin[k][expdof] = expbin[k];
                }
                for (int k=0; k<nomegas; k++)
                {
                    group_traj2[nn].omega[k][expdof] = expomega[k];
                    group_traj2[nn].bincounts[k][expdof] = expbincount[k];
                }
            }
        }
        for (int j=0; j<options.ndof; j++)
        {
            if (group_traj2[nn].coordinate.size() == 0)
            {
                std::cerr << "\nERROR: Could not find as many trajectories as degrees of freedom.  Exiting." << std::endl;
                std::exit(1);
            }
            else
            {
                ndof++;
            }
        }
        if (options.bVerbose)
        {
            std::cout << group_traj2[nn].nbins.size() << " DoF file" << grammar << " found for trajectory " << i << std::endl;
        }
        nn++;
    }
    group_traj = group_traj2;
    return;
}

void WHAM::cpp_wham_assign_bins()
{
    int this_bin;
    std::vector<int> bin_set(options.ndof);
    for (int i=0; i<nexp; i++)
    {
        int n = group_traj[i].experiment;
        if (options.bVerbose)
        {
            std::cout << "Doing bin assigments for experiment " << n << std::endl;
        }
        
        int ncoords = group_traj[i].coordinate.size();
        for (int j=0; j<ncoords; j++)
        {
            for (int k=0; k<options.ndof; k++)
            {
                /* This is SLOW */
                bin_set[k] = group_traj[i].frame_bin[j][k].bin1D;
            }
            this_bin = map1d(bin_set);
            for (int k=0; k<options.ndof; k++)
            {
                group_traj[i].frame_bin[j][k].binND = this_bin;
            }
        }
    }
    return;
}

void WHAM::cpp_wham_count_bins()
{
    int this_bin;
    for (int i=0; i<nexp; i++)
    {
        last_frame = group_traj[i].coordinate.size();
        if (last_frame > options.fN && options.fN != -1)
        {
            last_frame = options.fN;
        }
        if (last_frame < options.f0)
        {
            std::cerr << "\nError!  You chose to start analyzing on frame " << options.f0 << ", but there are only " << last_frame << " total frames in trajectory " << i << ".\n" << std::endl;
            std::exit(1);
        }
        int nf = 0;
        for (int j=options.f0; j<last_frame; j++)
        {
            this_bin = group_traj[i].frame_bin[j][0].binND;
            wham_args.counts[this_bin]++;
            nf += 1;
        }
        wham_args.sample[i] = nf;
    }

    if (options.bVerbose)
    {
        std::cout << "\n" << options.ndof << "-Dimensional State counts: [";
    }
    nzeros = 0;
    for (int i=0; i<nstates; i++)
    {
        if (options.bVerbose)
        {
            std::cout << wham_args.counts[i] << " ";
        }
        if (wham_args.counts[i] == 0)
        {
            nzeros++;
        }
    }
    if (options.bVerbose)
    {
        std::cout << "]" << std::endl;
        std::cout << "\nSamples: [ ";
        for (int i=0; i<nexp; i++)
        {
            std::cout << wham_args.sample[i] << " ";
        }
        std::cout << "]\n";
    }
    /* Report how many bins were unvisited */
    std::cout << "\n" << nzeros << "/" << nstates << " (" << nzeros*invnstates*100 << "%) bins unvisited.\n";
    return;
}

void WHAM::cpp_wham_make_ND_map()
{
    std::vector<t_map> m(nstates);
    std::vector<int> vND(options.ndof);
    map = m;
    for (int i=0; i<nstates; i++)
    {
        map[i].bin1D = vND;
        for (int j=0; j<options.ndof; j++)
        {
            if (j == options.ndof - 1)
            {
                vND[j]++;
            }
            if (vND[j] >= group_traj[0].nbins[j] && j > 0)
            {
                vND[j-1]++;
                vND[j] = 0;
            }
        }
    }
    if (options.bVerbose)
    {
        std::cout << "\nMapping " << options.ndof << "-dimensional array onto 1 dimensional array.\n1D <- D0 D1 D2..." << std::endl;
        for (int i=0 ; i<nstates; i++)
        {
            std::cout << i << " <- ";
            for (int j=0; j<options.ndof; j++)
            {
                std::cout << map[i].bin1D[j] << " ";
            }
            std::cout << std::endl;
        }
    }
    return;
}

void WHAM::cpp_wham_mapOmega()
{
    for (int h=0; h<nexp; h++)
    {
        std::vector<t_map> m(nstates);
        std::vector<int> vND(options.ndof);
        for (int i=0; i<nstates; i++)
        {
            for (int j=0; j<options.ndof; j++)
            {
                wham_args.omegas[h][i] *= group_traj[h].omega[vND[j]][j];
                if (j == options.ndof - 1)
                {
                    vND[j]++;
                }
                if (vND[j] >= group_traj[h].nbins[j] && j > 0)
                {
                    vND[j-1]++;
                    vND[j] = 0;
                }
            }
        }
    }

    if (options.bVerbose)
    {
        std::cout << "Built (" << nexp << "x" << nstates << ") omega matrix." << std::endl;
        for (int i=0; i<nexp; i++)
        {
            std::cout << "Omega for experiment " << i << ":\n [ ";
            for (int j=0; j<nstates; j++)
            {
                std::cout << wham_args.omegas[i][j] << " ";
            }
            std::cout << "]" << std::endl;
        }
    }
    for (int i=0; i<nexp; i+=1)
    {
        //wham_args.t[i] = group_traj[i][0].T();
        wham_args.t[i] = group_traj[i].t;
    }
    return;
}

void WHAM::cpp_wham_doWHAM()
{
    wham_args.nstates = nstates;
    opttrajectory.wham_init(wham_args, options);
    opt_prob = opttrajectory.PROB();
    opt_pmf = opttrajectory.PMF();
}

int WHAM::map1d(std::vector<int> bins)
{
    for (int i=0; i<nstates; i++)
    {
        bool match = false;
        for (int j=0; j<options.ndof; j++)
        {
            if (bins[j] == map[i].bin1D[j])
            {
                match = true;
            }
            else
            {
                match = false;
                break;
            }
        }
        if (match)
        {
            return i;
        }
    }
    std::cerr << "Unable to map ";
    int nbins = bins.size();
    for (int i=0; i<nbins; i++)
    {
        std::cerr << bins[i] << " ";
    }
    std::cerr << "to 1-Dimensional array.\n";
    std::exit(1);

    return -1;
}





int WHAM::NFrames()
{
    return nframes;
}
