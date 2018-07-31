#include "global_defs.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include <omp.h>

//#define prop_alpha_mean 0
//#define prop_alpha_var 5

using namespace std;

void jump_model_alpha()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_alpha;
    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_I /*reduction(+:nb_alpha_included)*/ private(r, A, old_alpha, new_phi,old_phi)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            old_alpha=alpha[i];

        // propose new alpha value
            if (!alpha_included[i])
                alpha[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_alpha[i],sqrt(var_alpha[i]));
            else
                alpha[i]=0;

        // change the state of alpha
            alpha_included[i]=!alpha_included[i];

        // calculate A
            A=0;

            for (int g=0;g<G;g++)
            {
                // calculate old and new value of theta
                new_phi=exp(-(alpha[i]+beta[g]));
                old_phi=exp(-(old_alpha+beta[g]));

                A+= gammaln(new_phi)-gammaln(old_phi)
                    - gammaln(new_phi*freq_ancestral[i]) + gammaln(old_phi*freq_ancestral[i])
                    - gammaln(new_phi*(1-freq_ancestral[i])) + gammaln(old_phi*(1-freq_ancestral[i]))
                    + freq_ancestral[i]*(new_phi-old_phi)*log(group[g].locus[i].p)
                    + (1-freq_ancestral[i])*(new_phi-old_phi)*log(1-group[g].locus[i].p);
            }

            if (alpha_included[i]) //if we add parameter
                A+= log_prior_alpha(alpha[i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
                    -(-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))-log(prior_odds);
        // inverse if we remove
            else
                A+= (-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))
                    -log_prior_alpha(alpha[i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

            r=randgen_parallel[omp_get_thread_num()].randDblExc();

        // reject proposed value
            if (log(r)>A)
            {
                alpha[i]=old_alpha;
                alpha_included[i]=!alpha_included[i];
            }
            /*else
            {
                if (alpha_included[i])
                    nb_alpha_included++;
                else
                    nb_alpha_included--;
            }*/
        }
    }
}


void jump_model_alpha_codominant()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_alpha;
    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_I /*reduction(+:nb_alpha_included)*/ private(r, A, old_alpha, new_phi,old_phi)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            old_alpha=alpha[i];

        // propose new alpha value
            if (!alpha_included[i])
                alpha[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_alpha[i],sqrt(var_alpha[i]));
            else
                alpha[i]=0;

        // change the state of alpha
            alpha_included[i]=!alpha_included[i];

        // calculate A
            A=0;

            double old_l=0;
            for (int g=0;g<G;g++)
            {
                // calculate old and new value of theta
                old_phi=exp(-(old_alpha+beta[g]));

                for (int k=0;k<group[g].locus[i].ar;k++)
                    old_l-=gammaln(old_phi*freq_locus[i].allele[k]);

                double temp=0;
                for (int k=0;k<group[g].locus[i].ar;k++)
                    temp+=old_phi*freq_locus[i].allele[k];

                old_l+=gammaln(temp);

                for (int k=0;k<group[g].locus[i].ar;k++)
                    old_l+=(old_phi*freq_locus[i].allele[k]-1)*log(group[g].locus[i].allele[k]);
            }

            double new_l=0;
            for (int g=0;g<G;g++)
            {
                // calculate old and new value of theta
                new_phi=exp(-(alpha[i]+beta[g]));

                for (int k=0;k<group[g].locus[i].ar;k++)
                    new_l-=gammaln(new_phi*freq_locus[i].allele[k]);

                double temp=0;
                for (int k=0;k<group[g].locus[i].ar;k++)
                    temp+=new_phi*freq_locus[i].allele[k];

                new_l+=gammaln(temp);

                for (int k=0;k<group[g].locus[i].ar;k++)
                    new_l+=(new_phi*freq_locus[i].allele[k]-1)*log(group[g].locus[i].allele[k]);
            }

            A=-old_l+new_l;

            if (alpha_included[i]) //if we add parameter
                A+= log_prior_alpha(alpha[i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
                    -(-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))-log(prior_odds);
        // inverse if we remove
            else
                A+= (-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))
                    -log_prior_alpha(alpha[i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

            r=randgen_parallel[omp_get_thread_num()].randDblExc();

        // reject proposed value
            if (log(r)>A)
            {
                alpha[i]=old_alpha;
                alpha_included[i]=!alpha_included[i];
            }
           /* else
            {
                if (alpha_included[i])
                    nb_alpha_included++;
                else
                    nb_alpha_included--;
            }*/
        }
    }
}

void jump_model_eta()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_eta;
    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_I private(r, A, old_eta, new_phi, old_phi)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int p=0;p<P;p++) // cycle over loci
            {

                old_eta=eta[p][i];

            // propose new alpha value
                if (!eta_included[p][i])
                    eta[p][i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_eta[p][i],sqrt(var_eta[p][i]));
                else
                    eta[p][i]=0;

            // change the state of alpha
                eta_included[p][i]=!eta_included[p][i];

            // calculate A
                A=0;
                for (int g=0;g<pressure[p].member.size();g++)
                {
                    int cur_group=pressure[p].member[g];

                    for (int j_g=0;j_g<group[cur_group].member.size();j_g++)
                    {
                        // calculate old and new value of theta
                        int cur_pop=group[cur_group].member[j_g];
                        new_phi=exp(-(eta[p][i]+theta[cur_pop]));
                        old_phi=exp(-(old_eta+theta[cur_pop]));

                        A+= gammaln(new_phi)-gammaln(old_phi)
                            - gammaln(new_phi*group[cur_group].locus[i].p) + gammaln(old_phi*group[cur_group].locus[i].p)
                            - gammaln(new_phi*(1-group[cur_group].locus[i].p)) + gammaln(old_phi*(1-group[cur_group].locus[i].p))
                            + group[cur_group].locus[i].p*(new_phi-old_phi)*log(pop[cur_pop].locus[i].p)
                            + (1-group[cur_group].locus[i].p)*(new_phi-old_phi)*log(1-pop[cur_pop].locus[i].p);
                    }
                }

                if (eta_included[p][i]) //if we add parameter
                    A+= log_prior_alpha(eta[p][i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
                        -(-0.5*log(2*M_PI*var_eta[p][i])-((eta[p][i]-mean_eta[p][i])*(eta[p][i]-mean_eta[p][i]))/(2*var_eta[p][i]))-log(prior_odds);
            // inverse if we remove
                else
                    A+= (-0.5*log(2*M_PI*var_eta[p][i])-((eta[p][i]-mean_eta[p][i])*(eta[p][i]-mean_eta[p][i]))/(2*var_eta[p][i]))
                        -log_prior_alpha(eta[p][i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

                r=randgen_parallel[omp_get_thread_num()].randDblExc();

            // reject proposed value
                if (log(r)>A)
                {
                    eta[p][i]=old_eta;
                    eta_included[p][i]=!eta_included[p][i];
                }
                /*else
                {
                    if (eta_included[g][i])
                        nb_eta_included[g]++;
                    else
                        nb_eta_included[g]--;
                }*/
            }
        }
    }
}
/*
void jump_model_eta_codominant()
{

    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_eta2;
    double new_phi,old_phi; // new and old values of phi

//    double old_log_likelihood; // old loglikelihood
    double diff_log_likelihood; // old loglikelihood

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, old_eta2, new_phi, old_phi, diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            //for (int p=0;p<P;p++) // cycle over loci
            int p=0;            
            p=randgen_parallel[omp_get_thread_num()].randInt(P-1);
            {
                old_eta2=eta2[p][i];

            // propose new alpha value
                if (!eta2_included[p][i])
                    eta2[p][i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_eta2[p][i],sqrt(var_eta2[p][i]));
                else
                    eta2[p][i]=0;

            // change the state of alpha
                eta2_included[p][i]=!eta2_included[p][i];

            // calculate A
                A=0;

                long double old_l=0;
                for (int g=0;g<pressure[p].member.size();g++)
                {
                    int cur_group=pressure[p].member[g];
                    for (int j_g=0;j_g<group[cur_group].member.size();j_g++)
                    {
                        int cur_pop=group[cur_group].member[j_g];
                        old_phi=exp(-(old_eta2+theta[cur_pop]));

                        old_l+=gammaln(old_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+old_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+old_phi*group[cur_group].locus[i].allele[k])
                                   -gammaln(old_phi*group[cur_group].locus[i].allele[k]);
                    }
                }

                long double new_l=0;
                for (int g=0;g<pressure[p].member.size();g++)
                {
                    int cur_group=pressure[p].member[g];
                    for (int j_g=0;j_g<group[cur_group].member.size();j_g++)
                    {
                        int cur_pop=group[cur_group].member[j_g];
                        new_phi=exp(-(eta2[p][i]+theta[cur_pop]));

                        new_l+=gammaln(new_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+new_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+new_phi*group[cur_group].locus[i].allele[k])
                                -gammaln(new_phi*group[cur_group].locus[i].allele[k]);
                    }
                }

            // store the old loglikelihood and calculate the new loglikelihood
                //old_log_likelihood=log_likelihood;
                diff_log_likelihood=-old_l+new_l;
                //log_likelihood=old_log_likelihood-old_l+new_l;

                A=diff_log_likelihood;

                if (eta2_included[p][i]) //if we add parameter
                    A+= log_prior_alpha(eta2[p][i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
                        -(-0.5*log(2*M_PI*var_eta2[p][i])-((eta2[p][i]-mean_eta2[p][i])*(eta2[p][i]-mean_eta2[p][i]))/(2*var_eta2[p][i]))-log(prior_odds);
            // inverse if we remove
                else
                    A+= (-0.5*log(2*M_PI*var_eta2[p][i])-((eta2[p][i]-mean_eta2[p][i])*(eta2[p][i]-mean_eta2[p][i]))/(2*var_eta2[p][i]))
                        -log_prior_alpha(eta2[p][i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

                r=randgen_parallel[omp_get_thread_num()].randDblExc();

            // reject proposed value
                if (log(r)>A)
                {
                    eta2[p][i]=old_eta2;
                    eta2_included[p][i]=!eta2_included[p][i];
                    //log_likelihood=old_log_likelihood;
                }
                else
                {
                    log_likelihood=log_likelihood+diff_log_likelihood;
                    //if (eta_included[g][i])
                    //    nb_eta_included[g]++;
                   // else
                   //     nb_eta_included[g]--;
                }
            }
        }
    }
}*/



void jump_model_eta_codominant_ok_old()
{

    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_eta2;
    double new_phi,old_phi; // new and old values of phi

//    double old_log_likelihood; // old loglikelihood
    double diff_log_likelihood; // old loglikelihood

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, old_eta2, new_phi, old_phi, diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int p=0;p<P;p++) 
            {
                
                int g=0;
                int group_index=0;
                bool all_present;
                do {
                    group_index=randgen_parallel[omp_get_thread_num()].randInt(pressure[p].member.size()-1);
                    g=pressure[p].member[group_index];
                    eta2_included[g][i]=!eta2_included[g][i];
                    all_present=true;
                    for (int g2=0;g2<pressure[p].member.size();g2++)
                        all_present=all_present && eta2_included[pressure[p].member[g2]][i];
                    eta2_included[g][i]=!eta2_included[g][i];
                } while (all_present && pressure[p].member.size()>1);

                old_eta2=eta2[g][i];

            // propose new alpha value
                if (!eta2_included[g][i])
                    eta2[g][i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_eta2[g][i],sqrt(var_eta2[g][i]));
                else
                    eta2[g][i]=0;

            // change the state of alpha
                eta2_included[g][i]=!eta2_included[g][i];

            // calculate A
                A=0;

                long double old_l=0;
                
                    for (int j_g=0;j_g<group[g].member.size();j_g++)
                    {
                        int cur_pop=group[g].member[j_g];
                        old_phi=exp(-(old_eta2+theta[cur_pop]));

                        old_l+=gammaln(old_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+old_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+old_phi*group[g].locus[i].allele[k])
                                   -gammaln(old_phi*group[g].locus[i].allele[k]);
                    }
                

                long double new_l=0;
               
                    for (int j_g=0;j_g<group[g].member.size();j_g++)
                    {
                        int cur_pop=group[g].member[j_g];
                        new_phi=exp(-(eta2[g][i]+theta[cur_pop]));

                        new_l+=gammaln(new_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+new_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+new_phi*group[g].locus[i].allele[k])
                                -gammaln(new_phi*group[g].locus[i].allele[k]);
                    }
                

            // store the old loglikelihood and calculate the new loglikelihood
                //old_log_likelihood=log_likelihood;
                diff_log_likelihood=-old_l+new_l;
                //log_likelihood=old_log_likelihood-old_l+new_l;

                A=diff_log_likelihood;

                if (eta2_included[g][i]) //if we add parameter
                    A+= log_prior_alpha(eta2[g][i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
                        -(-0.5*log(2*M_PI*var_eta2[g][i])-((eta2[g][i]-mean_eta2[g][i])*(eta2[g][i]-mean_eta2[g][i]))/(2*var_eta2[g][i]))-log(prior_odds);
            // inverse if we remove
                else
                    A+= (-0.5*log(2*M_PI*var_eta2[g][i])-((eta2[g][i]-mean_eta2[g][i])*(eta2[g][i]-mean_eta2[g][i]))/(2*var_eta2[g][i]))
                        -log_prior_alpha(eta2[g][i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

                r=randgen_parallel[omp_get_thread_num()].randDblExc();

            // reject proposed value
                if (log(r)>A)
                {
                    eta2[g][i]=old_eta2;
                    eta2_included[g][i]=!eta2_included[g][i];
                    //log_likelihood=old_log_likelihood;
                }
                else
                {
                    log_likelihood=log_likelihood+diff_log_likelihood;
                    /*if (eta_included[g][i])
                        nb_eta_included[g]++;
                    else
                        nb_eta_included[g]--;*/
                }
            }
        }
    }
}





void jump_model_eta_codominant()
{

    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_eta;
    double new_phi,old_phi; // new and old values of phi

//    double old_log_likelihood; // old loglikelihood
    double diff_log_likelihood; // old loglikelihood

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, old_eta, new_phi, old_phi, diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            // loop over pressures
            for (int p=0;p<P;p++) {
                // propose to move to  convergent evolution model in the pressure if
                //      - under convergent evolution model already
                // OR ( - more than one group in the pressure
                // AND  - all groups in the pressures are under selection in the current iteration
                // AND  - draw a random number > 0.5
                //     )
                // otherwise, do a normal move
                // random number to propose convergent evolution move
                double r2=randgen_parallel[omp_get_thread_num()].randDblExc();
                // check if all groups in the pressure are under selection
                int nb_present=0;
                for (int g=0;g<pressure[p].member.size();g++)
                    nb_present+=eta2_included[pressure[p].member[g]][i];

               /* bool all_present=true;
                for (int g=0;g<pressure[p].member.size();g++)                
                    all_present=all_present && eta2_included[pressure[p].member[g]][i];*/

                //if (i==150 && all_present) cout << r2 << endl;
                if (eta_included[p][i] || (pressure[p].member.size()>1 && nb_present==pressure[p].member.size()-1/*all_present*/ && r2>=0.5)) {
                    
                    /*f (i==150) {
                    cout << "propose CE model: ";
                    for (int g=0;g<pressure[p].member.size();g++)
                        cout << eta2_included[pressure[p].member[g]][i] << " ";
                    cout <<  eta_included[p][i];
                    cout << endl;
                    getchar();
                    }*/
                    // save the state
                    old_eta=eta[p][i];
                    
                    double* old_eta2=new double[pressure[p].member.size()];
                    for (int g=0;g<pressure[p].member.size();g++)
                        old_eta2[g]=eta2[pressure[p].member[g]][i];
                    
                    int g_null;

                    if (!eta_included[p][i]) { // propose to switch to convergent evolution model
                        // propose the new eta value as the average of eta2 values
                        eta[p][i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_eta[p][i],sqrt(var_eta[p][i]));
                       /* eta[p][i]=0;
                        for (int g=0;g<pressure[p].member.size();g++)
                            eta[p][i]+=eta2[pressure[p].member[g]][i];
                        eta[p][i]/=pressure[p].member.size()-1;*/
                        // set all eta2 to the same value to have correct calculations based on eta2 everywhere in the program without testing the current model
                        for (int g=0;g<pressure[p].member.size();g++)
                            eta2[pressure[p].member[g]][i]=eta[p][i];

                    }
                    else { // propose to switch to all groups under different selection coeff
                        eta[p][i]=0;
                        // propose the new coefficients for eta2 and change state
                        for (int g=0;g<pressure[p].member.size();g++)
                            eta2[pressure[p].member[g]][i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_eta2[pressure[p].member[g]][i],sqrt(var_eta2[pressure[p].member[g]][i]));
                    int group_index=0;
                    group_index=randgen_parallel[omp_get_thread_num()].randInt(pressure[p].member.size()-1);
                    g_null=pressure[p].member[group_index];
                    eta2[g_null][i]=0;
                    }

                    // change the state of eta and eta2
                    eta_included[p][i]=!eta_included[p][i];
                    for (int g=0;g<pressure[p].member.size();g++)
                        eta2_included[pressure[p].member[g]][i]=!eta2_included[pressure[p].member[g]][i];

                    // calculate A
                    A=0;
                    long double old_l=0;
                    long double new_l=0;
                    
                    for (int g_index=0;g_index<pressure[p].member.size();g_index++) {
                        int g=pressure[p].member[g_index];
                        
                        for (int j_g=0;j_g<group[g].member.size();j_g++)
                        {
                            int cur_pop=group[g].member[j_g];
                            if (!eta_included[p][i]) //if we switch to normal model
                                old_phi=exp(-(old_eta+theta[cur_pop]));
                            else old_phi=exp(-(old_eta2[g_index]+theta[cur_pop]));

                            old_l+=gammaln(old_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+old_phi);
                            for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                                old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+old_phi*group[g].locus[i].allele[k])
                                       -gammaln(old_phi*group[g].locus[i].allele[k]);
                        }

                        for (int j_g=0;j_g<group[g].member.size();j_g++)
                        {
                            int cur_pop=group[g].member[j_g];
                            new_phi=exp(-(eta2[g][i]+theta[cur_pop]));

                            new_l+=gammaln(new_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+new_phi);
                            for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                                new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+new_phi*group[g].locus[i].allele[k])
                                       -gammaln(new_phi*group[g].locus[i].allele[k]);
                        }
                    }

                    // store the old loglikelihood and calculate the new loglikelihood
                    //old_log_likelihood=log_likelihood;
                    diff_log_likelihood=-old_l+new_l;
                    //log_likelihood=old_log_likelihood-old_l+new_l;

                    A=diff_log_likelihood;
                    
                   /* if (eta_included[p][i]) {//if we switch to convergent evolution model
                        A+=log_prior_alpha(eta[p][i])+log(2)-log(pressure[p].member.size());
                        for (int g_index=0;g_index<pressure[p].member.size();g_index++) {
                            int g=pressure[p].member[g_index];
                            A+= ((int)(!eta2_included[g][i]))*(-0.5*log(2*M_PI*var_eta2[g][i])-((old_eta2[g_index]-mean_eta2[g][i])*(old_eta2[g_index]-mean_eta2[g][i]))/(2*var_eta2[g][i]))
                                -log_prior_alpha(old_eta2[g_index]);
                            A-= (-0.5*log(2*M_PI*var_eta[p][i])-((eta[p][i]-mean_eta[p][i])*(eta[p][i]-mean_eta[p][i]))/(2*var_eta[p][i]));

                        }
                    }
                    else {//if we switch to normal evolution model
                        A-=log_prior_alpha(old_eta)+log(2)-log(pressure[p].member.size());
                        for (int g_index=0;g_index<pressure[p].member.size();g_index++) {
                            int g=pressure[p].member[g_index];
                            A-= ((int)(g!=g_null))*(-0.5*log(2*M_PI*var_eta2[g][i])-((eta2[g][i]-mean_eta2[g][i])*(eta2[g][i]-mean_eta2[g][i]))/(2*var_eta2[g][i]))
                                -log_prior_alpha(eta2[g][i]);
                            A+= (-0.5*log(2*M_PI*var_eta[p][i])-((old_eta-mean_eta[p][i])*(old_eta-mean_eta[p][i]))/(2*var_eta[p][i]));
                        }
                    }*/


                    
                    if (eta_included[p][i]) {//if we switch to convergent evolution model
                        A+=log_prior_alpha(eta[p][i])+log(2)-log(pressure[p].member.size());
                        for (int g_index=0;g_index<pressure[p].member.size();g_index++) {
                            int g=pressure[p].member[g_index];
                            A+= ((int)(!eta2_included[g][i]))*(-0.5*log(2*M_PI*var_eta2[g][i])-((old_eta2[g_index]-mean_eta2[g][i])*(old_eta2[g_index]-mean_eta2[g][i]))/(2*var_eta2[g][i]))
                                    -((int)(!eta2_included[g][i]))*log_prior_alpha(old_eta2[g_index]);
                        }
                        A-= (-0.5*log(2*M_PI*var_eta[p][i])-((eta[p][i]-mean_eta[p][i])*(eta[p][i]-mean_eta[p][i]))/(2*var_eta[p][i]));                        
                    }
                    else {//if we switch to normal evolution model
                        A-=log_prior_alpha(old_eta)+log(2)-log(pressure[p].member.size());
                        for (int g_index=0;g_index<pressure[p].member.size();g_index++) {
                            int g=pressure[p].member[g_index];
                            A-= ((int)(g!=g_null))*(-0.5*log(2*M_PI*var_eta2[g][i])-((eta2[g][i]-mean_eta2[g][i])*(eta2[g][i]-mean_eta2[g][i]))/(2*var_eta2[g][i]))
                                    -((int)(g!=g_null))*log_prior_alpha(eta2[g][i]);
                        }
                        A+= (-0.5*log(2*M_PI*var_eta[p][i])-((old_eta-mean_eta[p][i])*(old_eta-mean_eta[p][i]))/(2*var_eta[p][i]));
                    }
                    
                    r=randgen_parallel[omp_get_thread_num()].randDblExc();

                // reject proposed value
                    if (log(r)>A)
                    {
                        eta[p][i]=old_eta;
                        for (int g=0;g<pressure[p].member.size();g++)
                            eta2[pressure[p].member[g]][i]=old_eta2[g];

                        eta_included[p][i]=!eta_included[p][i];
                        for (int g=0;g<pressure[p].member.size();g++)
                            eta2_included[pressure[p].member[g]][i]=!eta2_included[pressure[p].member[g]][i];
                        //log_likelihood=old_log_likelihood;
                    }
                    else
                    {
                       // if (i==150) cout << " -> accepted ";
                        log_likelihood=log_likelihood+diff_log_likelihood;
                        if (eta_included[p][i]) {
                            for (int g=0;g<pressure[p].member.size();g++)
                                eta2_included[pressure[p].member[g]][i]=false;
                        }
                        else {
                            eta2_included[g_null][i]=false;
                        }
                      /* if (i==150) {                     for (int g=0;g<pressure[p].member.size();g++)
                        cout << eta2_included[pressure[p].member[g]][i] << " ";
                    cout <<  eta_included[p][i];
                    cout << endl;}*/
                        //if (eta_included[g][i])
                        //    nb_eta_included[g]++;
                        //else
                        //    nb_eta_included[g]--;
                    }
                    delete[] old_eta2;
                } // end of convergent evolution move
                else {
                    // "normal" move: propose to switch one group in the pressure at random
                    int g=0;
                    int group_index=0;
                    bool all_present;
                    do {
                        group_index=randgen_parallel[omp_get_thread_num()].randInt(pressure[p].member.size()-1);
                        g=pressure[p].member[group_index];
                        eta2_included[g][i]=!eta2_included[g][i];
                        all_present=true;
                        for (int g2=0;g2<pressure[p].member.size();g2++)
                            all_present=all_present && eta2_included[pressure[p].member[g2]][i];
                        eta2_included[g][i]=!eta2_included[g][i];
                    } while (all_present && pressure[p].member.size()>1);                    
                    
                    /*if (i==150) {
                    cout << "propose normal move in group " << g << ": ";
                    for (int g=0;g<pressure[p].member.size();g++)
                        cout << eta2_included[pressure[p].member[g]][i] << " ";
                    cout <<  eta_included[p][i];
                    getchar();
                    }*/
                    double old_eta2=eta2[g][i];

                // propose new alpha value
                    if (!eta2_included[g][i])
                        eta2[g][i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_eta2[g][i],sqrt(var_eta2[g][i]));
                    else
                        eta2[g][i]=0;

                // change the state of alpha
                    eta2_included[g][i]=!eta2_included[g][i];

                // calculate A
                    A=0;

                    long double old_l=0;
                    for (int j_g=0;j_g<group[g].member.size();j_g++)
                    {
                        int cur_pop=group[g].member[j_g];
                        old_phi=exp(-(old_eta2+theta[cur_pop]));

                        old_l+=gammaln(old_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+old_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+old_phi*group[g].locus[i].allele[k])
                                   -gammaln(old_phi*group[g].locus[i].allele[k]);
                    }

                    long double new_l=0;
                    for (int j_g=0;j_g<group[g].member.size();j_g++)
                    {
                        int cur_pop=group[g].member[j_g];
                        new_phi=exp(-(eta2[g][i]+theta[cur_pop]));

                        new_l+=gammaln(new_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+new_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+new_phi*group[g].locus[i].allele[k])
                                   -gammaln(new_phi*group[g].locus[i].allele[k]);
                    }

                // store the old loglikelihood and calculate the new loglikelihood
                    //old_log_likelihood=log_likelihood;
                    diff_log_likelihood=-old_l+new_l;
                    //log_likelihood=old_log_likelihood-old_l+new_l;

                    A=diff_log_likelihood;

                    int nb_present=0;
                    for (int g=0;g<pressure[p].member.size();g++)
                        nb_present+=eta2_included[pressure[p].member[g]][i];

                    if (eta2_included[g][i]) {//if we add parameter
                        A+= log_prior_alpha(eta2[g][i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
                            -(-0.5*log(2*M_PI*var_eta2[g][i])-((eta2[g][i]-mean_eta2[g][i])*(eta2[g][i]-mean_eta2[g][i]))/(2*var_eta2[g][i]))-log(prior_odds);
                        if (pressure[p].member.size()>2 && (nb_present==pressure[p].member.size()-1))
                            A+=log(pressure[p].member.size())-log(pressure[p].member.size()-1);
                    }
                // inverse if we remove
                    else {
                        A+= (-0.5*log(2*M_PI*var_eta2[g][i])-((eta2[g][i]-mean_eta2[g][i])*(eta2[g][i]-mean_eta2[g][i]))/(2*var_eta2[g][i]))
                            -log_prior_alpha(eta2[g][i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));
                        if (pressure[p].member.size()>2 && (nb_present==pressure[p].member.size()-2))
                            A+=log(pressure[p].member.size()-1)-log(pressure[p].member.size());
                    }

                    r=randgen_parallel[omp_get_thread_num()].randDblExc();

                // reject proposed value
                    if (log(r)>A)
                    {
                        eta2[g][i]=old_eta2;
                        eta2_included[g][i]=!eta2_included[g][i];
                        //log_likelihood=old_log_likelihood;
                    }
                    else
                    {
                        //if (i==150) cout << " -> accepted ";
                        log_likelihood=log_likelihood+diff_log_likelihood;
                        //if (eta_included[g][i])
                        //    nb_eta_included[g]++;
                        //else
                        //    nb_eta_included[g]--;
                         /*if (i==150){                   for (int g=0;g<pressure[p].member.size();g++)
                        cout << eta2_included[pressure[p].member[g]][i] << " ";
                    cout <<  eta_included[p][i];
                    cout << endl;}*/
                    }
                } // end of normal move
            } // end of loop over pressures
        } // end of test if discarded locus
    } // end of loop over loci
} // end of function






