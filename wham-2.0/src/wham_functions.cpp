#include "wham_functions.h"

extern "C" {
    void sgemv_(const char *trans,
                const int *m,
                const int *n,
                const float *alpha,
                const float *a,
                const int *lda,
                const float *x,
                const int *incx,
                const float *beta,
                const float *y,
                const int *incy);
}
void sgemv(char trans,
           int m,
           int n,
           float alpha,
           float *a,
           int lda,
           float *x,
           int incx,
           float beta,
           float *y,
           int incy)
{
    return sgemv_(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}


void DOWHAM::wham_init(t_wham args, t_options options)
{
    wham_args = args;
    wham_options = options;
    nexperiments = options.ntraj / options.ndof;
    TransposeOmega();
    /* First WHAM step */
    /* WHAM iterations */
    std::vector<float> previous_step(wham_args.nstates,1.0/wham_args.nstates);
    std::vector<float> current_step(wham_args.nstates,1.0/wham_args.nstates);
    int points = 0;
    std::clock_t start = std::clock();
    while (points < wham_options.iter)
    {
        // slow step
        WhamStep(current_step);
        if ( square_diff(current_step,previous_step) < wham_options.tol)
        {
            std::cout << "Converged to " << wham_options.tol << " at step " << points << " after " << (std::clock() - start)/(float)(CLOCKS_PER_SEC/1000) << " ms." << std::endl;
            break;
        }
        /* 
         * Set the previous one as the result we just found and take a 
         * new step
         */
        else
        {
            previous_step = current_step;
        }
        points++;
    }
    if (wham_options.bVerbose)
    {
        std::cout << "\n[ ";
        for (int i=0; i<wham_args.nstates; i++)
        {
            std::cout << current_step[i] << " ";
        }
        std::cout << "]\n";
    }
    opt_trajectory = current_step;
    wham_pmf();
    wham_prob();
    return;
}

void DOWHAM::WhamStep(std::vector<float> &step)
{
#ifdef _OPENMP
    omp_set_num_threads(omp_get_max_threads());
#endif
    std::vector<float> fms(nexperiments,0);
    int vsize = wham_args.nstates; 
    std::vector<float> Momegas_times_trajectory(wham_args.nstates,0);
    sgemv('N',nexperiments,vsize,1,&c_major_omega[0],nexperiments,&step[0],1,1,&Momegas_times_trajectory[0],1);
    for (int i=0; i<nexperiments; i++)
    {
        fms[i] = wham_args.sample[i] / Momegas_times_trajectory[i];
    }
    vsize = nexperiments;
    std::vector<float> MomegasT_times_fms(wham_args.nstates,0);
    sgemv('N',wham_args.nstates,vsize,1,&c_major_omega_transpose[0],wham_args.nstates,&fms[0],1,1,&MomegasT_times_fms[0],1);
    for (int i=0; i<wham_args.nstates; i++)
    {
        step[i] = wham_args.counts[i] / MomegasT_times_fms[i];
    }
    vec_normalize(step);
    return;
}

void DOWHAM::TransposeOmega()
{
    std::vector<float> rows_by_columns(wham_args.omegas[0].size()*wham_args.omegas.size());
    c_major_omega = rows_by_columns;
    c_major_omega_transpose = rows_by_columns;
    int n = 0;
    for (int i=0; i<(int)wham_args.omegas[0].size(); i++)
    {
        for (int j=0; j<(int)wham_args.omegas.size(); j++)
        {
            c_major_omega[n] = (float)wham_args.omegas[j][i];       
            n++;
        }
    }
    
    n = 0;
    for (int i=0; i<(int)wham_args.omegas.size(); i++)
    {
        for (int j=0; j<(int)wham_args.omegas[0].size(); j++)
        {
            c_major_omega_transpose[n] = (float)wham_args.omegas[i][j];       
            n++;
        }
    }

    return;
}

float DOWHAM::square_diff(std::vector<float> a, std::vector<float> b)
{
    if (a.size() != b.size())
    {
        std::cerr << "\nERROR!  Size of vectors do not match, cannot find the square difference\n" << std::endl;
        std::exit(1);
    }
    std::vector<float> difference(a.size());
    for (int i=0; i<(int)a.size(); i++)
    {
        difference[i] = (a[i]-b[i])*(a[i]-b[i]);
    }
    return vec_sum(difference);
}

void DOWHAM::wham_pmf()
{
    /* Calculate the average temperature */
    float avg_experiment_T, sum;
    for (int i=0; i<(float)wham_args.sample.size();i++)
    {
        avg_experiment_T += wham_args.sample[i] * wham_args.t[i];
        sum += wham_args.sample[i];
    }
    avg_experiment_T /= sum;
    /* PMF_bin(i) = -kT*ln(w(i)).  Assume that any bin that had 0 visists had a pmf of 1000 */
    std::vector<float> p(wham_args.nstates);
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
    std::vector<float> p(wham_args.nstates);
    probability = p;
    float beta = 1./(KB*wham_options.t);
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

std::vector<float> DOWHAM::PMF()
{
    return potential;
}

std::vector<float> DOWHAM::PROB()
{
    return probability;
}
