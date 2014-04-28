#include "wham_functions.h"

void DOWHAM::wham_init(t_wham args, t_options options)
{
    wham_args = args;
    wham_options = options;
    nexperiments = options.ntraj / options.ndof;
    TransposeOmega();
    std::vector<double> w0(wham_args.nstates,1.0/wham_args.nstates);
    
    /* First WHAM step */
    std::vector<std::vector<double> > t;
    std::vector<double> ptraj_init(wham_args.nstates);
    trajectory = t;
    trajectory.push_back(w0);
    /* WHAM iterations */
    int points = 0;
    while (points < wham_options.iter)
    {
        trajectory.push_back(ptraj_init);
        trajectory[points+1] = WhamStep(points);
        if (square_diff(trajectory[points+1],trajectory[points]) < wham_options.tol)
        {
            std::cout << "Converged to " << wham_options.tol << " at step " << points << std::endl;
            break;
        }
        points++;
    }
    if (wham_options.bVerbose)
    {
        std::cout << "\n[ ";
        for (int i=0; i<wham_args.nstates; i++)
        {
            std::cout << trajectory[points][i] << " ";
        }
        std::cout << "]\n";
    }
    opt_trajectory = trajectory[points];
    wham_pmf();
    wham_prob();
    return;
}

std::vector<double> DOWHAM::WhamStep(int point)
{
    std::vector<double> fms(nexperiments),wguess(wham_args.nstates);
    for (int i=0; i<nexperiments; i++)
    {
        fms[i] = wham_args.sample[i]/vec_dot(trajectory[point],wham_args.omegas[i]);
    }
    for (int i=0; i<wham_args.nstates; i++)
    {
        wguess[i] = wham_args.counts[i]/vec_dot(omegaT[i],fms);
    }
    vec_normalize(wguess);
    return wguess;
}

void DOWHAM::TransposeOmega()
{
    std::vector<double> ncol(wham_args.omegas.size());
    std::vector<std::vector<double> > oT(wham_args.omegas[0].size());
    omegaT = oT;
    for (int i=0; i<wham_args.omegas[0].size(); i++)
    {
        omegaT[i] = ncol;
        for (int j=0; j<wham_args.omegas.size();j++)
        {
            omegaT[i][j] = wham_args.omegas[j][i];
        }
    }
    return;
}

double DOWHAM::square_diff(std::vector<double> a, std::vector<double> b)
{
    if (a.size() != b.size())
    {
        std::cerr << "\nERROR!  Size of vectors do not match, cannot find the square difference\n" << std::endl;
        exit(1);
    }
    std::vector<double> difference(a.size());
    for (int i=0; i<a.size(); i++)
    {
        difference[i] = (a[i]-b[i])*(a[i]-b[i]);
    }
    return vec_sum(difference);
}

void DOWHAM::wham_pmf()
{
    /* Calculate the average temperature */
    double avg_experiment_T = vec_dot(wham_args.sample, wham_args.t)/vec_sum(wham_args.sample);
    /* PMF_bin(i) = -kT*ln(w(i)).  Assume that any bin that had 0 visists had a pmf of 1000 */
    std::vector<double> p(wham_args.nstates);
    potential = p;
    for (int i=0; i<wham_args.nstates; i++)
    {
        
        if (opt_trajectory[i] <= 1e-30)
        {
            potential[i] = 1000;
        }
        else
        {
            potential[i] = -KB * avg_experiment_T * log(opt_trajectory[i]);
        }
    }
    return;
}

void DOWHAM::wham_prob()
{
    std::vector<double> p(wham_args.nstates);
    probability = p;
    double beta = 1./(KB*wham_options.t);
    for (int i=0; i<wham_args.nstates; i++)
    {
        probability[i] = exp(-beta*potential[i]);
    }
    if (wham_options.bVerbose && true == false)
    {
        for (int i=0;i<wham_args.nstates;i++)
        {
            std::cout << i << " " << probability[i] << std::endl;
        }
    }
}

std::vector<double> DOWHAM::PMF()
{
    return potential;
}

std::vector<double> DOWHAM::PROB()
{
    return probability;
}
