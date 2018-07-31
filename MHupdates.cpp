#include "global_defs.h"
#include <iostream>
#include <algorithm>
#include <math.h>

#include <omp.h>

#define sd_prior_beta 1 // 1.8 sd of the normal prior for beta
#define mean_prior_beta -1 // -2 mean of the normal prior for beta

#define sd_prior_theta 1 // 1.8 sd of the normal prior for beta
#define mean_prior_theta -1 // -2 mean of the normal prior for beta

#define lambda 1 // parameter for the non informative dirichlet prior of allele frequencies

#define gamma 10 // parameter for the variance between delta and mu

#define epsilon 1e-6 // limit for allele frequencies (epsilon,1-epsilon)

#define PI 3.14159265358979

using namespace std;

//////////////////////////////////
// update ancestral allele frequencies
///////////////////////////////////
void update_ancestral_freq()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // forward and backward proposition probability

    double old_freq; // store old allele freq if the move is rejected

    double u; // proposal random value

    //omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for SCHED_I private(r, A, FwPrp, BwPrp, old_freq, u)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            // store old value
            old_freq=freq_ancestral[i];

            // calculate forward proposition probability =q(p',p)
            BwPrp=min(1.0-epsilon,old_freq+e_ancestral[i])-max(0.0+epsilon,old_freq-e_ancestral[i]);

            // propose the new values
            u=randgen_parallel[omp_get_thread_num()].randDblExc(BwPrp)+max(0.0+epsilon,old_freq-e_ancestral[i]);
            if (u<=epsilon) u=epsilon;
            if (u>=1-epsilon) u=1-epsilon;
            freq_ancestral[i]=u;

            // calculate backward proposition probability =q(p,p')
            FwPrp=min(1.0-epsilon,freq_ancestral[i]+e_ancestral[i])-max(0.0+epsilon,freq_ancestral[i]-e_ancestral[i]);

            // ratio A
            A= log(BwPrp) - log(FwPrp);
            /*+ (a_p-1)*log((1-freq_ancestral[i])/(1-old_freq))+(a_p-1)*log(freq_ancestral[i]/old_freq);*/
            double phi;
            for (int g=0;g<G;g++)
            {
                phi=exp(-(alpha[i]+beta[g]));
                A+= gammaln(phi*old_freq) + gammaln(phi*(1-old_freq))
                    - gammaln(phi*freq_ancestral[i]) - gammaln(phi*(1-freq_ancestral[i]))
                    + phi*(freq_ancestral[i]-old_freq)*(log(group[g].locus[i].p) - log(1-group[g].locus[i].p));
            }
            
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                freq_ancestral[i]=old_freq;// restore old value of allele frequency
            }
            else
            {
                acc_freq_ancestral[i]++;
            }
        }
    } // end of cycle over loci

}


/////////////////////////////////////
//  update a_p (beta prior for p)
/////////////////////////////////////
void update_a_p()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double u;
    double e;

    double old_a; // old value of a

    double m=0;
    double s=1;

    // keep the old value
    old_a=a_p;

    // make the move
    e=randgen.randNorm(0,var_prop_a_p);
    a_p=a_p*exp(e);

    A= 0.5*((log(old_a)-m)/s)*((log(old_a)-m)/s)
       - 0.5*((log(a_p)-m)/s)*((log(a_p)-m)/s);

    A+=I*(gammaln(2.0*a_p)+2.0*gammaln(old_a)-gammaln(2.0*old_a)-2.0*gammaln(a_p));
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i])
            A+=(a_p-old_a)*(log(freq_ancestral[i])+log(1-freq_ancestral[i]));
    }


    r=randgen.randDblExc();
    // reject proposed value
    if (log(r)>A)
    {
        a_p=old_a; // restore old value of sigma_square
    }
    else
    {
        acc_a_p++;
    }

}


//////////////////////////////////
// update allele frequencies
///////////////////////////////////
void update_freq()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // forward and backward proposition probability
    double FwPrp2,BwPrp2; // forward and backward proposition probability

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_freq; // store old allele freq if the move is rejected

    //double log_likelihood2;

    double old_g,new_g; // store new and old phenotype frequency

    double u; // proposal random value

    double g;
    
    #pragma omp parallel for SCHED_J reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, FwPrp2, BwPrp2, diff_log_likelihood, old_freq, old_g, new_g, u, g)
    for (int j=0;j<J;j++) //cycle over populations
    {
        for (int i=0;i<I;i++) // cycle over loci
        {
            if (!discarded_loci[i])
            {
                double e_small=e_freq[i][j]/100;
                double prop_small=0.5;

                // store old value
                old_freq=pop[j].locus[i].p;

                old_g= pop[j].locus[i].p*pop[j].locus[i].p
                       +2*pop[j].locus[i].p*(1-pop[j].locus[i].p)*(1-f[j])
                       +f[j]*pop[j].locus[i].p*(1-pop[j].locus[i].p);

                // calculate forward proposition probability =q(p',p)
                BwPrp=min(1.0-epsilon,old_freq+e_freq[i][j])-max(0.0+epsilon,old_freq-e_freq[i][j]);
                BwPrp2=min(1.0-epsilon,old_freq+e_small)-max(0.0+epsilon,old_freq-e_small);
                // propose the new values
                    if (randgen_parallel[omp_get_thread_num()].rand()<1-prop_small)
                        u=randgen_parallel[omp_get_thread_num()].rand(BwPrp)+max(0.0+epsilon,old_freq-e_freq[i][j]);
                    else
                        u=randgen_parallel[omp_get_thread_num()].rand(BwPrp2)+max(0.0+epsilon,old_freq-e_small);

                if (u<=epsilon) u=epsilon;
                if (u>=1-epsilon) u=1-epsilon;
                pop[j].locus[i].p=u;

                // calculate backward proposition probability =q(p,p')
                FwPrp=min(1.0-epsilon,pop[j].locus[i].p+e_freq[i][j])-max(0.0+epsilon,pop[j].locus[i].p-e_freq[i][j]);
                FwPrp2=min(1.0-epsilon,pop[j].locus[i].p+e_small)-max(0.0+epsilon,pop[j].locus[i].p-e_small);
                // store the old logikelihood and calculate the new logikelihood
                //old_log_likelihood=log_likelihood;

                new_g= pop[j].locus[i].p*pop[j].locus[i].p
                       +2*pop[j].locus[i].p*(1-pop[j].locus[i].p)*(1-f[j])
                       +f[j]*pop[j].locus[i].p*(1-pop[j].locus[i].p);

                diff_log_likelihood= pop[j].locus[i].nA1*(log(new_g)-log(old_g))
                                     + (pop[j].locus[i].n-pop[j].locus[i].nA1)*(log(1-new_g)-log(1-old_g));

                //#pragma omp critical(likelihood)
                //log_likelihood=old_log_likelihood+diff_log_likelihood;

                double phi_sc=exp(-(eta2[pop[j].group][i]+theta[j]));
                // ratio A
                A= diff_log_likelihood
                   +(phi_sc*group[pop[j].group].locus[i].p-1)*(log(pop[j].locus[i].p)-log(old_freq))
                   + (phi_sc*(1-group[pop[j].group].locus[i].p)-1)*(log(1-pop[j].locus[i].p)-log(1-old_freq))
                   + log((1-prop_small)*BwPrp+prop_small*BwPrp2) - log((1-prop_small)*FwPrp+prop_small*FwPrp2);
                
                r=randgen_parallel[omp_get_thread_num()].randDblExc();
                // reject proposed value
                if (log(r)>A)
                {
                    //log_likelihood=old_log_likelihood;  // restore likelihood
                    pop[j].locus[i].p=old_freq;// restore old value of allele frequency
                }
                else
                {
                    acc_freq[i][j]++;                
                	log_likelihood=log_likelihood+diff_log_likelihood;
                }
            }
        } // end of cycle over loci


    } // end of cycle over populations

}



//////////////////////////////////
// update group allele frequencies
///////////////////////////////////
void update_freq_group()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // forward and backward proposition probability

//double old_log_likelihood; // old logikelihood
//double diff_log_likelihood; // difference of the logikelihoods
    double old_freq; // store old allele freq if the move is rejected

//double log_likelihood2;

    double u; // proposal random value

//double theta; // store current theta value

    //omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for SCHED_I  private(r, A, FwPrp, BwPrp, old_freq, u)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int g=0;g<G;g++) //cycle over populations
            {
                // store old value
                old_freq=group[g].locus[i].p;

                // calculate forward proposition probability =q(p',p)
                BwPrp=min(1.0-epsilon,old_freq+e_freq_group[i][g])-max(0.0+epsilon,old_freq-e_freq_group[i][g]);

                // propose the new values
                u=randgen_parallel[omp_get_thread_num()].rand(BwPrp)+max(0.0+epsilon,old_freq-e_freq_group[i][g]);
                if (u<=epsilon) u=epsilon;
                if (u>=1-epsilon) u=epsilon;
                group[g].locus[i].p=u;

                // calculate backward proposition probability =q(p,p')
                FwPrp=min(1.0-epsilon,group[g].locus[i].p+e_freq_group[i][g])-max(0.0+epsilon,group[g].locus[i].p-e_freq_group[i][g]);

                A=0;
                double phi_sc;
                for (int j=0;j<group[g].member.size();j++)
                {
                    int cur_pop=group[g].member[j];
                    phi_sc=exp(-(eta2[g][i]+theta[cur_pop]));
                    A+= gammaln(phi_sc*old_freq) + gammaln(phi_sc*(1-old_freq))
                        - gammaln(phi_sc*group[g].locus[i].p) - gammaln(phi_sc*(1-group[g].locus[i].p))
                        + phi_sc*(group[g].locus[i].p-old_freq)*(log(pop[cur_pop].locus[i].p) - log(1-pop[cur_pop].locus[i].p));
                }

                // ratio A
                double phi_ct=exp(-(alpha[i]+beta[g]));
                A+= (phi_ct*freq_ancestral[i]-1)*(log(group[g].locus[i].p)-log(old_freq))
                    + (phi_ct*(1-freq_ancestral[i])-1)*(log(1-group[g].locus[i].p)-log(1-old_freq))
                    + log(BwPrp) - log(FwPrp);

                r=randgen_parallel[omp_get_thread_num()].randDblExc();
                // reject proposed value
                if (log(r)>A)
                {
//                  log_likelihood=old_log_likelihood;  // restore likelihood
                    group[g].locus[i].p=old_freq;// restore old value of allele frequency
                }
                else
                {
                    acc_freq_group[i][g]++;
                }
            }
        } // end of cycle over loci


    } // end of cycle over populations

}


//////////////////////////////////
// update fis randomly (NOT ESTIMATED)
///////////////////////////////////
void update_f_random()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q
    double old_f; // store old f if the move is rejected
    double u; // proposal random value
    if (prior_fis_unif)
    {
        //omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for SCHED_J private(r, A, FwPrp, BwPrp, u,old_f)
        for (int j=0;j<J;j++) // cycle over populations
        {
            // store old value
            old_f=f[j];
            // calculate forward proposition probability =q(p',p)
            BwPrp=min(prior_fis_hb,old_f+e_f[j])-max(prior_fis_lb,old_f-e_f[j]);
            // propose the new values
            u=randgen_parallel[omp_get_thread_num()].randDblExc(BwPrp)+max(prior_fis_lb,old_f-e_f[j]);
            f[j]=u;
            // calculate backward proposition probability =q(p,p')
            FwPrp=min(1.0,f[j]+e_f[j])-max(0.0,f[j]-e_f[j]);
            // ratio A
            A= log(BwPrp) - log(FwPrp);
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)	 f[j]=old_f;// restore old value of allele frequency
            else
            {
                //log_likelihood=allelecount_loglikelihood();
                acc_f[j]++;
            }
        }
        log_likelihood=allelecount_loglikelihood();
    }
    else   // beta prior
    {
        #pragma omp parallel for  private(r, A, FwPrp, BwPrp, u,old_f)
        for (int j=0;j<J;j++) // cycle over populations
        {
            // store old value
            old_f=f[j];
            // calculate forward proposition probability =q(p',p)
            BwPrp=min(1.0,old_f+e_f[j])-max(0.0,old_f-e_f[j]);
            // propose the new values
            u=randgen_parallel[omp_get_thread_num()].rand(BwPrp)+max(0.0,old_f-e_f[j]);
            if (u<=0.0001) u=0.0001;
            if (u>=0.9999) u=0.9999;
            f[j]=u;
            // calculate backward proposition probability =q(p,p')
            FwPrp=min(1.0,f[j]+e_f[j])-max(0.0,f[j]-e_f[j]);
            A= (prior_fis_a-1)*(log(f[j])-log(old_f))
               + (prior_fis_b-1)*(log(1-f[j])-log(1-old_f))
               + log(BwPrp) - log(FwPrp);
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)	 f[j]=old_f;// restore old value of allele frequency
            else
            {
                //log_likelihood=allelecount_loglikelihood();
                acc_f[j]++;
            }
        }
        log_likelihood=allelecount_loglikelihood();
    }
}

//////////////////////////////////
// update theta
///////////////////////////////////

void update_theta()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double BwPrp,FwPrp;
    double u;

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of theta

    double old_theta;

    //omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for SCHED_J private(r, A, FwPrp, BwPrp, e, u,new_phi,old_phi,old_theta)
    for (int j=0;j<J;j++) // cycle over populations
    {
        // store old value
        old_theta=theta[j];
        // make the move
        e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_theta[j]);
        theta[j]=theta[j]+e;

        // ratio A
        A=(old_theta-theta[j])*(old_theta+theta[j]-2*mean_prior_theta)/(2*sd_prior_theta*sd_prior_theta);

        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                // calculate old and new value of theta
                new_phi=exp(-(eta2[pop[j].group][i]+theta[j]));
                old_phi=exp(-(eta2[pop[j].group][i]+old_theta));
                A+= gammaln(new_phi)-gammaln(old_phi)
                    - gammaln(new_phi*group[pop[j].group].locus[i].p) + gammaln(old_phi*group[pop[j].group].locus[i].p)
                    - gammaln(new_phi*(1-group[pop[j].group].locus[i].p)) + gammaln(old_phi*(1-group[pop[j].group].locus[i].p))
                    + group[pop[j].group].locus[i].p*(new_phi-old_phi)*log(pop[j].locus[i].p)
                    + (1-group[pop[j].group].locus[i].p)*(new_phi-old_phi)*log(1-pop[j].locus[i].p);
            }
        }
		
        r=randgen_parallel[omp_get_thread_num()].randDblExc();
        // reject proposed value
        if (log(r)>A)
        {
            theta[j]=old_theta;// restore old value of beta
        }
        else
        {
            acc_theta[j]++;
        }

    } // end of cycle over populations

}

//////////////////////////////////
// update beta
///////////////////////////////////
void update_beta()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_beta; // store old allele freq if the move is rejected
    double BwPrp,FwPrp;
    double u;

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of theta

//    #pragma omp parallel for SCHED_G private(r, A, FwPrp, BwPrp, e, u,new_phi,old_phi,old_beta)
    for (int g=0;g<G;g++) // cycle over populations
    {
        // store old value
        old_beta=beta[g];

        // make the move
        e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_beta[g]);
        beta[g]=beta[g]+e;

        // ratio A
        A=(old_beta-beta[g])*(old_beta+beta[g]-2*mean_prior_beta)/(2*sd_prior_beta*sd_prior_beta);

       #pragma omp parallel for SCHED_I reduction(+:A) private(new_phi,old_phi)
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                // calculate old and new value of theta
                new_phi=exp(-(alpha[i]+beta[g]));
                old_phi=exp(-(alpha[i]+old_beta));

                A+= gammaln(new_phi)-gammaln(old_phi)
                    - gammaln(new_phi*freq_ancestral[i]) + gammaln(old_phi*freq_ancestral[i])
                    - gammaln(new_phi*(1-freq_ancestral[i])) + gammaln(old_phi*(1-freq_ancestral[i]))
                    + freq_ancestral[i]*(new_phi-old_phi)*log(group[g].locus[i].p)
                    + (1-freq_ancestral[i])*(new_phi-old_phi)*log(1-group[g].locus[i].p);
            }
        }

        r=randgen_parallel[omp_get_thread_num()].randDblExc();
        // reject proposed value
        if (log(r)>A)
        {
            beta[g]=old_beta;// restore old value of beta
        }
        else
        {
            acc_beta[g]++;
        }

    } // end of cycle over populations

}

//////////////////////////////////
// update alpha i
///////////////////////////////////
void update_alpha()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_alpha; // store old allele freq if the move is rejected

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of theta
    
    #pragma omp parallel for SCHED_I /*reduction(+:alpha_updates)*/ private(r, A, old_alpha, e,new_phi,old_phi)
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i] && alpha_included[i])
        {
            //alpha_updates=alpha_updates+1;
            // store old value
            old_alpha=alpha[i];

            // make the move
            //std::cout << omp_get_thread_num();
            e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_alpha[i]);
            alpha[i]=alpha[i]+e;

            // ratio A
            A= (old_alpha*old_alpha-alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha);
        //        A=log_prior_alpha(alpha[i])-log_prior_alpha(old_alpha);
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

            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                alpha[i]=old_alpha;// restore old value of alpha
            }
            else
            {
                acc_alpha[i]++;
            }
        }
    }

}


//////////////////////////////////
// update eta i g
///////////////////////////////////
void update_eta()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_eta; // store old eta if the move is rejected

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of theta

    #pragma omp parallel for SCHED_I private(r, A, old_eta, e, new_phi, old_phi)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int p=0;p<P;p++) //cycle over populations
            {
                if (eta_included[p][i])
                {
                    //eta_updates[g]++;
                    // store old value
                    old_eta=eta[p][i];

                    // make the move
                    e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_eta[p][i]);
                    eta[p][i]=eta[p][i]+e;

                    // ratio A
                    A= (old_eta*old_eta-eta[p][i]*eta[p][i])/(2*sd_prior_alpha*sd_prior_alpha);

                    for (int g=0;g<pressure[p].member.size();g++)
                    {
                        int cur_group=pressure[p].member[g];
                        for (int j=0;j<group[cur_group].member.size();j++)
                        {
                        // calculate old and new value of theta
                            int cur_pop=group[cur_group].member[j];
                            new_phi=exp(-(eta[p][i]+theta[cur_pop]));
                            old_phi=exp(-(old_eta+theta[cur_pop]));

                            A+= gammaln(new_phi)-gammaln(old_phi)
                                - gammaln(new_phi*group[cur_group].locus[i].p) + gammaln(old_phi*group[cur_group].locus[i].p)
                                - gammaln(new_phi*(1-group[cur_group].locus[i].p)) + gammaln(old_phi*(1-group[cur_group].locus[i].p))
                               + group[cur_group].locus[i].p*(new_phi-old_phi)*log(pop[cur_pop].locus[i].p)
                               + (1-group[cur_group].locus[i].p)*(new_phi-old_phi)*log(1-pop[cur_pop].locus[i].p);
                        }
                    }

                    r=randgen_parallel[omp_get_thread_num()].randDblExc();
                    // reject proposed value
                    if (log(r)>A)
                    {
                        eta[p][i]=old_eta;// restore old value of alpha
                    }
                    else
                    {
                        acc_eta[p][i]++;
                    }
                }
            }
        }
    }
}



//////////////////////////////////
// update eta2 i g
///////////////////////////////////
void update_eta2()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_eta2; // store old eta if the move is rejected

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of theta

    #pragma omp parallel for SCHED_I private(r, A, old_eta2, e, new_phi, old_phi)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int g=0;g<G;g++) //cycle over populations
            {
                if (eta2_included[g][i])
                {
                    //eta_updates[g]++;
                    // store old value
                    old_eta2=eta2[g][i];

                    // make the move
                    e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_eta2[g][i]);
                    eta2[g][i]=eta2[g][i]+e;

                    // ratio A
                    A= (old_eta2*old_eta2-eta2[g][i]*eta2[g][i])/(2*sd_prior_alpha*sd_prior_alpha);

                    for (int j=0;j<group[g].member.size();j++)
                    {
                        // calculate old and new value of theta
                        int cur_pop=group[g].member[j];
                        new_phi=exp(-(eta2[g][i]+theta[cur_pop]));
                        old_phi=exp(-(old_eta2+theta[cur_pop]));

                        A+= gammaln(new_phi)-gammaln(old_phi)
                            - gammaln(new_phi*group[g].locus[i].p) + gammaln(old_phi*group[g].locus[i].p)
                            - gammaln(new_phi*(1-group[g].locus[i].p)) + gammaln(old_phi*(1-group[g].locus[i].p))
                            + group[g].locus[i].p*(new_phi-old_phi)*log(pop[cur_pop].locus[i].p)
                            + (1-group[g].locus[i].p)*(new_phi-old_phi)*log(1-pop[cur_pop].locus[i].p);
                    }

                    r=randgen_parallel[omp_get_thread_num()].randDblExc();
                    // reject proposed value
                    if (log(r)>A)
                    {
                        eta2[g][i]=old_eta2;// restore old value of alpha
                    }
                    else
                    {
                        acc_eta2[g][i]++;
                    }
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////
// For codominant markers
/////////////////////////////////////////////


//////////////////////////////////
// update allele frequencies in groups
///////////////////////////////////
void update_freq_group_codominant()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

    double old_m,old_n; // store old allele freq if the move is rejected

    int m,n; // alleles to be changed, chosen at random

    //double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihood
//double log_likelihood2;

    double u; // proposal random value

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, old_m,old_n, m,n,u,diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int g=0;g<G;g++) // cycle over groups
            {
                if (freq_locus[i].ar==2) //case of SNPs, no need to choose allele at random
                {
                    m=0;
                    n=1;
                }
                else
                {
                    // choose at random the two allele frequences to change
                    m=randgen_parallel[omp_get_thread_num()].randInt(freq_locus[i].ar-1);
                    do
                    {
                        n=randgen_parallel[omp_get_thread_num()].randInt(freq_locus[i].ar-1);
                    }
                    while (n==m);
                }

                // store old values
                old_m=group[g].locus[i].allele[m];
                old_n=group[g].locus[i].allele[n];

                // store old modify part of likelihood
                // old_changed locus
                double old_l=0;
                for (int j=0;j<group[g].member.size();j++)
                {
                    int cur_pop=group[g].member[j];
                    double phi_sc=exp(-(eta2[g][i]+theta[cur_pop]));
                    old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[m]+phi_sc*group[g].locus[i].allele[m])
                           -factln(pop[cur_pop].locus[i].data_allele_count[m])-gammaln(phi_sc*group[g].locus[i].allele[m]);
                    old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[n]+phi_sc*group[g].locus[i].allele[n])
                           -factln(pop[cur_pop].locus[i].data_allele_count[n])-gammaln(phi_sc*group[g].locus[i].allele[n]);

                }
                // calculate inverse of forward and backward proposition probability =1/q
                BwPrp=min(old_m+old_n,old_m+e_freq_group[i][g])-max(0.0,old_m-e_freq_group[i][g]);

                // propose the new values
                u=randgen_parallel[omp_get_thread_num()].randDblExc(BwPrp)+max(0.0,old_m-e_freq_group[i][g]);
                group[g].locus[i].allele[m]=u;
                group[g].locus[i].allele[n]=old_n+(old_m-group[g].locus[i].allele[m]);


                FwPrp=min(group[g].locus[i].allele[m]+group[g].locus[i].allele[n],group[g].locus[i].allele[m]+e_freq_group[i][g])-max(0.0,group[g].locus[i].allele[m]-e_freq_group[i][g]);


                // store the old logikelihood and calculate the new logikelihood
                //old_log_likelihood=log_likelihood;

                double new_l=0;
                // calculate new changed locus
                for (int j=0;j<group[g].member.size();j++)
                {
                    int cur_pop=group[g].member[j];
                    double phi_sc=exp(-(eta2[g][i]+theta[cur_pop]));
                    new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[m]+phi_sc*group[g].locus[i].allele[m])
                           -factln(pop[cur_pop].locus[i].data_allele_count[m])-gammaln(phi_sc*group[g].locus[i].allele[m]);
                    new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[n]+phi_sc*group[g].locus[i].allele[n])
                           -factln(pop[cur_pop].locus[i].data_allele_count[n])-gammaln(phi_sc*group[g].locus[i].allele[n]);
                }
                //log_likelihood=old_log_likelihood-old_l+new_l;
                diff_log_likelihood=-old_l+new_l;
                // ratio A
                double phi_ct=exp(-(alpha[i]+beta[g]));
                A= -old_l+new_l + log(BwPrp) - log(FwPrp)
                   // bug found : was lambda*... -> (lambda-1)*... is correct
                   //	+ (lambda-1.0)*(log(group[g].locus[i].allele[m])+log(group[g].locus[i].allele[n])-log(old_m)-log(old_n));
                   + (phi_ct*freq_locus[i].allele[m]-1.0)* ( log(group[g].locus[i].allele[m])-log(old_m) )
                   + (phi_ct*freq_locus[i].allele[n]-1.0)* ( log(group[g].locus[i].allele[n])-log(old_n) );
                r=randgen_parallel[omp_get_thread_num()].randDblExc();
                // reject proposed value
                if (log(r)>A)
                {
                    //log_likelihood=old_log_likelihood;  // restore likelihood
                    group[g].locus[i].allele[m]=old_m;  // and old values of allele frequencies
                    group[g].locus[i].allele[n]=old_n;
                }
                else
                {
                    acc_freq_group[i][g]++;
                    log_likelihood=log_likelihood+diff_log_likelihood;
                }
            }
        } // end of cycle over loci
    } // end of cycle over groups

}



//////////////////////////////////
// update ancetral allele frequencies
///////////////////////////////////
void update_ancestral_freq_codominant()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

    double old_m,old_n; // store old allele freq if the move is rejected

    int m,n; // alleles to be changed, chosen at random

    double phi;

    double u; // proposal random value

    #pragma omp parallel for SCHED_I private(r, A, FwPrp, BwPrp, old_m,old_n, m,n,phi,u)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            if (freq_locus[i].ar==2) //case of SNPs, no need to choose allele at random
            {
                m=0;
                n=1;
            }
            else
            {
                // choose at random the two allele frequences to change
                m=randgen_parallel[omp_get_thread_num()].randInt(freq_locus[i].ar-1);
                do
                {
                    n=randgen_parallel[omp_get_thread_num()].randInt(freq_locus[i].ar-1);
                }
                while (n==m);
            }

            // store old values
            old_m=freq_locus[i].allele[m];
            old_n=freq_locus[i].allele[n];

            // store old modify part of likelihood
            // old_changed locus
            double old_l=0;
            for (int g=0;g<G;g++)
            {
                phi=exp(-(alpha[i]+beta[g]));
                for (int k=0;k<group[g].locus[i].ar;k++)
                    old_l-=gammaln(phi*freq_locus[i].allele[k]);

                double temp=0;
                for (int k=0;k<group[g].locus[i].ar;k++)
                    temp+=phi*freq_locus[i].allele[k];

                old_l+=gammaln(temp);

                for (int k=0;k<group[g].locus[i].ar;k++)
                    old_l+=(phi*freq_locus[i].allele[k]-1)*log(group[g].locus[i].allele[k]);
            }

            // calculate inverse of forward and backward proposition probability =1/q
            BwPrp=min(old_m+old_n,old_m+e_ancestral[i])-max(0.0,old_m-e_ancestral[i]);

            // propose the new values
            u=randgen_parallel[omp_get_thread_num()].randDblExc(BwPrp)+max(0.0,old_m-e_ancestral[i]);
            freq_locus[i].allele[m]=u;
            freq_locus[i].allele[n]=old_n+(old_m-freq_locus[i].allele[m]);

            FwPrp=min(freq_locus[i].allele[m]+freq_locus[i].allele[n],freq_locus[i].allele[m]+e_ancestral[i])-max(0.0,freq_locus[i].allele[m]-e_ancestral[i]);

            double new_l=0;
            // calculate new changed locus
            for (int g=0;g<G;g++)
            {
                phi=exp(-(alpha[i]+beta[g]));
                for (int k=0;k<group[g].locus[i].ar;k++)
                    new_l-=gammaln(phi*freq_locus[i].allele[k]);

                double temp=0;
                for (int k=0;k<group[g].locus[i].ar;k++)
                    temp+=phi*freq_locus[i].allele[k];

                new_l+=gammaln(temp);

                for (int k=0;k<group[g].locus[i].ar;k++)
                    new_l+=(phi*freq_locus[i].allele[k]-1)*log(group[g].locus[i].allele[k]);
            }

            // ratio A
            A= -old_l+new_l + log(BwPrp) - log(FwPrp)
               // bug found : was lambda*... -> (lambda-1)*... is correct
               + (lambda-1.0)*(log(freq_locus[i].allele[m])+log(freq_locus[i].allele[n])-log(old_m)-log(old_n));

            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                freq_locus[i].allele[m]=old_m;  // and old values of allele frequencies
                freq_locus[i].allele[n]=old_n;
            }
            else
            {
                acc_freq_ancestral[i]++;
            }
        }
    } // end of cycle over loci

}


//////////////////////////////////
// update alpha i
///////////////////////////////////
void update_alpha_codominant()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_alpha; // store old allele freq if the move is rejected

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_I /*reduction(+:alpha_updates)*/ private(r, A, old_alpha, e,new_phi,old_phi)
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i] && alpha_included[i])
        {
            //alpha_updates=alpha_updates+1;
            // store old value
            old_alpha=alpha[i];

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

            // make the move
            e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_alpha[i]);
            alpha[i]=alpha[i]+e;

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

            // ratio A
            A= -old_l+new_l +(old_alpha*old_alpha-alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha);

            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                alpha[i]=old_alpha;// restore old value of alpha
            }
            else
            {
                acc_alpha[i]++;
            }
        }
    }

}



//////////////////////////////////
// update eta i g
///////////////////////////////////
void update_eta_codominant()
{

    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_eta; // store old eta if the move is rejected

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood;

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, old_eta, e, new_phi, old_phi, diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int p=0;p<P;p++) //cycle over populations
            {
                if (eta_included[p][i])
                {                    
                    //eta_updates[g]++;
                    // store old value                    
                    old_eta=eta[p][i];

                    double old_l=0;
                    for (int g=0;g<pressure[p].member.size();g++)
                    {
                        int cur_group=pressure[p].member[g];
                        for (int j=0;j<group[cur_group].member.size();j++)
                        {
                           int cur_pop=group[cur_group].member[j];
                            old_phi=exp(-(old_eta+theta[cur_pop]));

                            old_l+=gammaln(old_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+old_phi);
                            for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                                old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+old_phi*group[cur_group].locus[i].allele[k])
                                       -gammaln(old_phi*group[cur_group].locus[i].allele[k]);
                        }
                    }

                    // make the move
                    e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_eta[p][i]);
                    eta[p][i]=eta[p][i]+e;

                    double new_l=0;
                    for (int g=0;g<pressure[p].member.size();g++)
                    {
                        int cur_group=pressure[p].member[g];
                        for (int j=0;j<group[cur_group].member.size();j++)
                        {
                            int cur_pop=group[cur_group].member[j];
                            new_phi=exp(-(eta[p][i]+theta[cur_pop]));

                            new_l+=gammaln(new_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+new_phi);
                            for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                                new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+new_phi*group[cur_group].locus[i].allele[k])
                                       -gammaln(new_phi*group[cur_group].locus[i].allele[k]);
                        }
                    }

                    // store the old logikelihood and calculate the new logikelihood
                    //old_log_likelihood=log_likelihood;

                    //log_likelihood=old_log_likelihood-old_l+new_l;
                    diff_log_likelihood=-old_l+new_l;
                    // ratio A
                    A= diff_log_likelihood +(old_eta*old_eta-eta[p][i]*eta[p][i])/(2*sd_prior_alpha*sd_prior_alpha);

                    r=randgen_parallel[omp_get_thread_num()].randDblExc();
                   
                    // reject proposed value
                    if (log(r)>A)
                    {
                        eta[p][i]=old_eta;// restore old value of alpha
                        //log_likelihood=old_log_likelihood;
                    }
                    else
                    {
                        log_likelihood=log_likelihood+diff_log_likelihood;
                        for (int g=0;g<pressure[p].member.size();g++)
                            eta2[pressure[p].member[g]][i]=eta[p][i];
                        acc_eta[p][i]++;
                    }
                    
                }
            }
        }
    }
}


//////////////////////////////////
// update eta2 i g
///////////////////////////////////
void update_eta2_codominant()
{

    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_eta2; // store old eta if the move is rejected

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood;

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, old_eta2, e, new_phi, old_phi, diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            for (int g=0;g<G;g++) //cycle over populations
            {
                if (eta2_included[g][i])
                {
                    //eta_updates[g]++;
                    // store old value
                    old_eta2=eta2[g][i];

                    double old_l=0;
                    for (int j=0;j<group[g].member.size();j++)
                    {
                        int cur_pop=group[g].member[j];
                        old_phi=exp(-(old_eta2+theta[cur_pop]));

                        old_l+=gammaln(old_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+old_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            old_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+old_phi*group[g].locus[i].allele[k])
                                   -gammaln(old_phi*group[g].locus[i].allele[k]);
                    }

                    // make the move
                    e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_eta2[g][i]);
                    eta2[g][i]=eta2[g][i]+e;

                    double new_l=0;
                    for (int j=0;j<group[g].member.size();j++)
                    {
                        int cur_pop=group[g].member[j];
                        new_phi=exp(-(eta2[g][i]+theta[cur_pop]));

                        new_l+=gammaln(new_phi)-gammaln(pop[cur_pop].locus[i].alleleCount+new_phi);
                        for (int k=0;k<pop[cur_pop].locus[i].ar;k++)
                            new_l+=gammaln(pop[cur_pop].locus[i].data_allele_count[k]+new_phi*group[g].locus[i].allele[k])
                                   -gammaln(new_phi*group[g].locus[i].allele[k]);
                    }

                    // store the old logikelihood and calculate the new logikelihood
                    //old_log_likelihood=log_likelihood;

                    //log_likelihood=old_log_likelihood-old_l+new_l;
                    diff_log_likelihood=-old_l+new_l;
                    // ratio A
                    A= diff_log_likelihood +(old_eta2*old_eta2-eta2[g][i]*eta2[g][i])/(2*sd_prior_alpha*sd_prior_alpha);

                    r=randgen_parallel[omp_get_thread_num()].randDblExc();
                    // reject proposed value
                    if (log(r)>A)
                    {
                        eta2[g][i]=old_eta2;// restore old value of alpha
                        //log_likelihood=old_log_likelihood;
                    }
                    else
                    {
                        log_likelihood=log_likelihood+diff_log_likelihood;
                        acc_eta2[g][i]++;
                    }
                }
            }
        }
    }

}


//////////////////////////////////
// update beta
///////////////////////////////////
void update_beta_codominant()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

    double old_beta; // store old allele freq if the move is rejected

    double e; // proposal random value

    double new_phi,old_phi; // new and old values of phi

//   #pragma omp parallel for SCHED_G private(r, A, e, new_phi,old_phi,old_beta)
    for (int g=0;g<G;g++)
    {
        // store old value
        old_beta=beta[g];

        double old_l=0;
        #pragma omp parallel for SCHED_I reduction(+:old_l) private(old_phi)
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                // calculate old and new value of theta
                old_phi=exp(-(alpha[i]+old_beta));

                for (int k=0;k<group[g].locus[i].ar;k++)
                    old_l+=-gammaln(old_phi*freq_locus[i].allele[k]);

                double temp=0;
                for (int k=0;k<group[g].locus[i].ar;k++)
                    temp+=old_phi*freq_locus[i].allele[k];

                old_l+=gammaln(temp);

                for (int k=0;k<group[g].locus[i].ar;k++)
                    old_l+=(old_phi*freq_locus[i].allele[k]-1)*log(group[g].locus[i].allele[k]);
            }
        }

        // make the move
        e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_beta[g]);
        beta[g]=beta[g]+e;

        double new_l=0;
        #pragma omp parallel for SCHED_I reduction(+:new_l) private(new_phi)
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                // calculate old and new value of theta
                new_phi=exp(-(alpha[i]+beta[g]));

                for (int k=0;k<group[g].locus[i].ar;k++)
                    new_l+=-gammaln(new_phi*freq_locus[i].allele[k]);

                double temp=0;
                for (int k=0;k<group[g].locus[i].ar;k++)
                    temp+=new_phi*freq_locus[i].allele[k];

                new_l+=gammaln(temp);

                for (int k=0;k<group[g].locus[i].ar;k++)
                    new_l+=(new_phi*freq_locus[i].allele[k]-1)*log(group[g].locus[i].allele[k]);
            }
        }

        // ratio A
        A= -old_l+new_l +
           (old_beta-beta[g])*(old_beta+beta[g]-2*mean_prior_beta)/(2*sd_prior_beta*sd_prior_beta);

        r=randgen_parallel[omp_get_thread_num()].randDblExc();
        // reject proposed value
        if (log(r)>A)
        {
            beta[g]=old_beta;// restore old value of beta
        }
        else
        {
            acc_beta[g]++;
        }

    } // end of cycle over populations

}



//////////////////////////////////
// update theta
///////////////////////////////////
void update_theta_codominant()
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)

//    double old_log_likelihood;
    double diff_log_likelihood;

    double e; // proposal random value

    double old_theta; // new and old values of theta

    double new_phi,old_phi; // new and old values of phi

    #pragma omp parallel for SCHED_J reduction(+:log_likelihood)  private(r, A, e, new_phi,old_phi,old_theta,diff_log_likelihood)
    for (int j=0;j<J;j++)
    {

        // store old value
        old_theta=theta[j];

        double old_l=0;
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                old_phi=exp(-(eta2[pop[j].group][i]+old_theta));
                old_l+=gammaln(old_phi)-gammaln(pop[j].locus[i].alleleCount+old_phi);
                for (int k=0;k<pop[j].locus[i].ar;k++)
                    old_l+=gammaln(pop[j].locus[i].data_allele_count[k]+old_phi*group[pop[j].group].locus[i].allele[k])
                           -gammaln(old_phi*group[pop[j].group].locus[i].allele[k]);
            }
        }

        // make the move
        e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_theta[j]);
        theta[j]=theta[j]+e;

        double new_l=0;
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
               new_phi=exp(-(eta2[pop[j].group][i]+theta[j]));
                new_l+=gammaln(new_phi)-gammaln(pop[j].locus[i].alleleCount+new_phi);
                for (int k=0;k<pop[j].locus[i].ar;k++)
                    new_l+=gammaln(pop[j].locus[i].data_allele_count[k]+new_phi*group[pop[j].group].locus[i].allele[k])
                           -gammaln(new_phi*group[pop[j].group].locus[i].allele[k]);
            }
        }

        // store the old logikelihood and calculate the new logikelihood
        //old_log_likelihood=log_likelihood;

        //log_likelihood=old_log_likelihood-old_l+new_l;
        diff_log_likelihood=-old_l+new_l;

        // ratio A
        A= diff_log_likelihood +
           (old_theta-theta[j])*(old_theta+theta[j]-2*mean_prior_theta)/(2*sd_prior_theta*sd_prior_theta);;

        r=randgen_parallel[omp_get_thread_num()].randDblExc();
        // reject proposed value
        if (log(r)>A)
        {
            theta[j]=old_theta;// restore old value of beta
            //log_likelihood=old_log_likelihood;
        }
        else
        {
            log_likelihood=log_likelihood+diff_log_likelihood;
            acc_theta[j]++;
        }

    } // end of cycle over populations

}


////////////////////////////////
// update f for intensity data
////////////////////////////////
void update_f_intensity() // selection2
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_f; // store old f if the move is rejected

    double old_log_likelihood_part; // old part logikelihood modified
    double new_log_likelihood_part; // new part logikelihood modified

    double u; // proposal random value

    #pragma omp parallel for SCHED_J reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, u,old_f,diff_log_likelihood,old_log_likelihood_part,new_log_likelihood_part)
    for (int j=0;j<J;j++)
    {

        // store old value
        old_f=f[j];

        // calculate forward proposition probability =q(p',p)
        BwPrp=min(1.0,old_f+e_f[j])-max(0.0,old_f-e_f[j]);

        // calculate part of the likelihood that will be modified
        old_log_likelihood_part=0;
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    old_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            }
        }

        // propose the new values
        u=randgen_parallel[omp_get_thread_num()].randDblExc(BwPrp)+max(0.0,old_f-e_f[j]);
        f[j]=u;

        // calculate backward proposition probability =q(p,p')
        FwPrp=min(1.0,f[j]+e_f[j])-max(0.0,f[j]-e_f[j]);

        // calculate part of the likelihood modified
        new_log_likelihood_part=0;
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    new_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            }
        }

        // store the old logikelihood and calculate the new logikelihood
        //old_log_likelihood=log_likelihood;

        diff_log_likelihood = new_log_likelihood_part-old_log_likelihood_part;

        //log_likelihood=old_log_likelihood+diff_log_likelihood;

        // ratio A
        A= diff_log_likelihood + log(BwPrp) - log(FwPrp);

        r=randgen_parallel[omp_get_thread_num()].randDblExc();
        // reject proposed value
        if (log(r)>A)
        {
            //log_likelihood=old_log_likelihood;  // restore likelihood
            f[j]=old_f;// restore old value of allele frequency
        }
        else
        {
            log_likelihood=log_likelihood+diff_log_likelihood;
            acc_f[j]++;
        }
    }
}

void update_freq_intensity()// selection2
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // forward and backward proposition probability

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_freq; // store old allele freq if the move is rejected

    double old_log_likelihood_part; // old part logikelihood modified
    double new_log_likelihood_part; // new part logikelihood modified

    double u; // proposal random value

    double g;

    #pragma omp parallel for SCHED_J reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, u,old_freq,diff_log_likelihood,old_log_likelihood_part,new_log_likelihood_part,g)
    for (int j=0;j<J;j++) //cycle over populations
    {
        for (int i=0;i<I;i++) // cycle over loci
        {
            if (!discarded_loci[i])
            {

                // store old value
                old_freq=pop[j].locus[i].p;

                // calculate forward proposition probability =q(p',p)
                BwPrp=min(1.0-epsilon,old_freq+e_freq[i][j])-max(0.0+epsilon,old_freq-e_freq[i][j]);

                // calculate part of the likelihood that will be modified
                old_log_likelihood_part=0;
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    old_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);

                // propose the new values
                u=randgen_parallel[omp_get_thread_num()].rand(BwPrp)+max(0.0+epsilon,old_freq-e_freq[i][j]);
                if (u<=epsilon) u=epsilon;
                if (u>=1-epsilon) u=epsilon;
                pop[j].locus[i].p=u;

                // calculate backward proposition probability =q(p,p')
                FwPrp=min(1.0-epsilon,pop[j].locus[i].p+e_freq[i][j])-max(0.0+epsilon,pop[j].locus[i].p-e_freq[i][j]);

                // calculate new modified part of the likelihood
                new_log_likelihood_part=0;
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    new_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);

                // store the old logikelihood and calculate the new logikelihood
                //old_log_likelihood=log_likelihood;

                diff_log_likelihood= new_log_likelihood_part-old_log_likelihood_part;

                //log_likelihood=old_log_likelihood+diff_log_likelihood;

                /*log_likelihood2=allelecount_logikelihood(Application);
                if (log_likelihood!=log_likelihood2)
                log_likelihood=log_likelihood2;
                */

                // ratio A
                double phi_sc=exp(-(eta[pop[j].group][i]+theta[j]));
                A= diff_log_likelihood
                   +(phi_sc*group[pop[j].group].locus[i].p-1)*(log(pop[j].locus[i].p)-log(old_freq))
                   + (phi_sc*(1-group[pop[j].group].locus[i].p)-1)*(log(1-pop[j].locus[i].p)-log(1-old_freq))
                   + log(BwPrp) - log(FwPrp);

                // maybe important bug : in following part theta and freq_ancestral[i]
                // are for non hierarchical version only !!! (forgot to change this in the first hierarchical version)
                //theta=exp(-(alpha[i]+beta[j]));
                // ratio A
                //A= diff_log_likelihood
                //   +(theta*freq_ancestral[i]-1)*(log(pop[j].locus[i].p)-log(old_freq))
                //   + (theta*(1-group[pop[j].group].locus[i].p)-1)*(log(1-pop[j].locus[i].p)-log(1-old_freq))
                //   + log(BwPrp) - log(FwPrp);

                r=randgen_parallel[omp_get_thread_num()].randDblExc();
                // reject proposed value
                if (log(r)>A)
                {
                    //log_likelihood=old_log_likelihood;  // restore likelihood
                    pop[j].locus[i].p=old_freq;// restore old value of allele frequency
                }
                else
                {
                    log_likelihood=log_likelihood+diff_log_likelihood;
                    acc_freq[i][j]++;
                }
            }
        } // end of cycle over loci


    } // end of cycle over populations

}

void update_mu_intensity()// selection2
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_mu; // store old mu if the move is rejected

    double old_log_likelihood_part; // old part logikelihood modified
    double new_log_likelihood_part; // new part logikelihood modified

    double u; // proposal random value
    //double s; // choose type of proposal

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, u,old_mu,diff_log_likelihood,old_log_likelihood_part,new_log_likelihood_part)
    for (int i = 0; i < I; i++)
    {
    
        if (!discarded_loci[i])
        {
            // store old value
            old_mu=mu[i];

            // calculate forward proposition probability =q(p',p)
//	BwPrp=min(1.0,old_mu+e_mu)-max(0.0,old_mu-e_mu);
            BwPrp=atan((1-old_mu)/e_mu[i]) - atan(abscence_pc-old_mu/e_mu[i]);

            // calculate part of the likelihood that will be modified
            old_log_likelihood_part=locus_likelihood[i];//0;
            /*for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    old_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            } */

            // propose the new values
//	u=randgen.randDblExc(BwPrp)+max(0.0,old_mu-e_mu);
//	mu[i]=u;
            do
            {
                mu[i]=old_mu;
                u=randgen_parallel[omp_get_thread_num()].randDblExc();
                mu[i]=mu[i]+e_mu[i]*tan(PI*(u-0.5));
            }
            while (mu[i]<=abscence_pc || mu[i]>=1);

            FwPrp=atan((1-mu[i])/e_mu[i]) - atan(abscence_pc-mu[i]/e_mu[i]);

            // calculate backward proposition probability =q(p,p')
//	FwPrp=min(1.0,mu[i]+e_mu)-max(0.0,mu[i]-e_mu);

            // calculate new modified part of the likelihood
            new_log_likelihood_part=0;
            for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    new_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            }

            // store the old logikelihood and calculate the new logikelihood
            //old_log_likelihood=log_likelihood;

            diff_log_likelihood=new_log_likelihood_part-old_log_likelihood_part;

            //log_likelihood=old_log_likelihood+diff_log_likelihood;


            // ratio A  to improve
            // cauchy :
            //A= diff_log_likelihood + log((delta[i]-old_mu)*(delta[i]-old_mu)+gamma*gamma) -  log((delta[i]-mu[i])*(delta[i]-mu[i])+gamma*gamma) + log(BwPrp) - log(FwPrp);
            // normal :
            A= diff_log_likelihood + gamma*gamma*(delta[i]-old_mu)*(delta[i]-old_mu)/(2*old_mu*old_mu) - gamma*gamma*(delta[i]-mu[i])*(delta[i]-mu[i])/(2*mu[i]*mu[i])  + log(BwPrp) - log(FwPrp);
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                //log_likelihood=old_log_likelihood;  // restore likelihood
                mu[i]=old_mu;// restore old value of allele frequency
            }
            else
            {
                log_likelihood=log_likelihood+diff_log_likelihood;
                locus_likelihood[i]=new_log_likelihood_part;
                acc_mu[i]++;
            }
        }
    } // end of loop over loci

}


void update_delta_intensity()// selection2
{
    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_delta; // store old delta if the move is rejected

    double old_log_likelihood_part; // old part logikelihood modified
    double new_log_likelihood_part; // new part logikelihood modified

    double u; // proposal random value
//    double s; // choose type of proposal

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, u,old_delta,diff_log_likelihood,old_log_likelihood_part,new_log_likelihood_part)
    for (int i = 0; i < I; i++)
    {   
        if (!discarded_loci[i])
        {
            // store old value
            old_delta=delta[i];

            // calculate forward proposition probability =q(p',p)
            //BwPrp=min(1.0,old_delta+e_delta)-max(0.0,old_delta-e_delta);
            BwPrp=atan((1-old_delta)/e_delta[i]) - atan(-old_delta/e_delta[i]);

            // calculate part of the likelihood that will be modified
            old_log_likelihood_part=locus_likelihood[i];
            /*for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    old_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            } */

            // propose the new values
//	u=randgen.randDblExc(BwPrp)+max(0.0,old_delta-e_delta);
//	delta[i]=u;

            do
            {
                delta[i]=old_delta;
                u=randgen_parallel[omp_get_thread_num()].randDblExc();
                delta[i]=delta[i]+e_delta[i]*tan(PI*(u-0.5));
            }
            while (delta[i]<=0 || delta[i]>=1);

            FwPrp=atan((1-delta[i])/e_delta[i]) - atan(-delta[i]/e_delta[i]);

            // calculate backward proposition probability =q(p,p')
            //FwPrp=min(1.0,delta[i]+e_delta)-max(0.0,delta[i]-e_delta);

            // calculate new modified part of the likelihood
            new_log_likelihood_part=0;
            for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    new_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            }

            // store the old logikelihood and calculate the new logikelihood
            //old_log_likelihood=log_likelihood;

            diff_log_likelihood=new_log_likelihood_part-old_log_likelihood_part;

            //log_likelihood=old_log_likelihood+diff_log_likelihood;


            // ratio A  to improve
            //cauchy:
            //A= diff_log_likelihood + log((old_delta-mu[i])*(old_delta-mu[i])+gamma*gamma) -  log((delta[i]-mu[i])*(delta[i]-mu[i])+gamma*gamma) + log(BwPrp)- log(FwPrp);
            // normal:
            A= diff_log_likelihood + gamma*gamma*(old_delta-delta[i])*(old_delta+delta[i]-2*mu[i])/(2*mu[i]*mu[i]) + log(BwPrp)- log(FwPrp);
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                //log_likelihood=old_log_likelihood;  // restore likelihood
                delta[i]=old_delta;// restore old value of allele frequency
            }
            else
            {
                log_likelihood=log_likelihood+diff_log_likelihood;
                locus_likelihood[i]=new_log_likelihood_part;
                acc_delta[i]++;
            }
        }
    } // end of loop over loci

}


void update_sigma1_intensity()// selection2
{

    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

//    double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_sigma1; // store old sigma1 if the move is rejected

    double old_log_likelihood_part; // old part logikelihood modified
    double new_log_likelihood_part; // new part logikelihood modified

//double u; // proposal random value
    double e; // for proposal step

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, e,old_sigma1,diff_log_likelihood,old_log_likelihood_part,new_log_likelihood_part)
    for (int i = 0; i < I; i++)
    {
        if (!discarded_loci[i])
        {

            // store old value
            old_sigma1=sigma1[i];

            // calculate forward proposition probability =q(p',p)
            //BwPrp=min(1.0,old_sigma1+var_prop_sigma1[i])-max(0.0,old_sigma1-var_prop_sigma1[i]);

            // calculate part of the likelihood that will be modified
            old_log_likelihood_part=locus_likelihood[i];
            /*for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    old_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            } */

            // propose the new values
            //u=randgen.randDblExc(BwPrp)+max(0.0,old_sigma1-var_prop_sigma1[i]);
            //sigma1[i]=u;

            // calculate backward proposition probability =q(p,p')
            //FwPrp=min(1.0,sigma1[i]+var_prop_sigma1[i])-max(0.0,sigma1[i]-var_prop_sigma1[i]);


            // make the move
            e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_sigma1[i]);
            sigma1[i]=sigma1[i]*exp(e);
			if (sigma1[i]<=epsilon) sigma1[i]=epsilon;
			
            // calculate new modified part of the likelihood
            new_log_likelihood_part=0;
            for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    new_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            }

            // store the old logikelihood and calculate the new logikelihood
            //old_log_likelihood=log_likelihood;
//	log_likelihood=allelecount_logikelihood(Application);

            diff_log_likelihood= new_log_likelihood_part-old_log_likelihood_part;

            //log_likelihood=old_log_likelihood+diff_log_likelihood;

            // ratio A  to improve
//	A= diff_log_likelihood + log(BwPrp) - log(FwPrp);
            A= diff_log_likelihood + log(sigma1[i]) - log(old_sigma1) + (old_sigma1 - sigma1[i])/0.05;
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                //log_likelihood=old_log_likelihood;  // restore likelihood
                sigma1[i]=old_sigma1;// restore old value of allele frequency
            }
            else
            {
                log_likelihood=log_likelihood+diff_log_likelihood;
                locus_likelihood[i]=new_log_likelihood_part;
                acc_sigma1[i]++;
            }
        }
    } // end of loop over loci

}


void update_sigma2_intensity()// selection2
{

    double r; // random value to accept/reject the move
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double FwPrp,BwPrp; // inverse of forward and backward proposition probability =1/q

    //double old_log_likelihood; // old logikelihood
    double diff_log_likelihood; // difference of the logikelihoods
    double old_sigma2; // store old sigma2 if the move is rejected

    double old_log_likelihood_part; // old part logikelihood modified
    double new_log_likelihood_part; // new part logikelihood modified

//double u; // proposal random value
    double e; // for proposal step

    #pragma omp parallel for SCHED_I reduction(+:log_likelihood) private(r, A, FwPrp, BwPrp, e,old_sigma2,diff_log_likelihood,old_log_likelihood_part,new_log_likelihood_part)
    for (int i = 0; i < I; i++)
    {
        if (!discarded_loci[i])
        {

            // store old value
            old_sigma2=sigma2[i];

            // calculate forward proposition probability =q(p',p)
            //BwPrp=min(1.0,old_sigma2+var_prop_sigma2[i])-max(0.0,old_sigma2-var_prop_sigma2[i]);

            // calculate part of the likelihood that will be modified
            old_log_likelihood_part=locus_likelihood[i];
            /*for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    old_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            } */

            // propose the new values
            //u=randgen.randDblExc(BwPrp)+max(0.0,old_sigma2-var_prop_sigma2[i]);
            //sigma2[i]=u;

            // calculate backward proposition probability =q(p,p')
            //FwPrp=min(1.0,sigma2[i]+var_prop_sigma2[i])-max(0.0,sigma2[i]-var_prop_sigma2[i]);


            // make the move
            e=randgen_parallel[omp_get_thread_num()].randNorm(0,var_prop_sigma2[i]);
            sigma2[i]=sigma2[i]*exp(e);
			if (sigma2[i]<=epsilon) sigma2[i]=epsilon;
			
            // calculate new modified part of the likelihood
            new_log_likelihood_part=0;
            for (int j = 0; j < J; j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                    new_log_likelihood_part+=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
            }

            // store the old logikelihood and calculate the new logikelihood
            //old_log_likelihood=log_likelihood;
//	log_likelihood=allelecount_logikelihood(Application);

            diff_log_likelihood= new_log_likelihood_part-old_log_likelihood_part;

            //log_likelihood=old_log_likelihood+diff_log_likelihood;

            // ratio A  to improve
//	A= diff_log_likelihood + log(BwPrp) - log(FwPrp);
            A= diff_log_likelihood + log(sigma2[i]) - log(old_sigma2) + (old_sigma2 - sigma2[i])/0.05;
            r=randgen_parallel[omp_get_thread_num()].randDblExc();
            // reject proposed value
            if (log(r)>A)
            {
                //log_likelihood=old_log_likelihood;  // restore likelihood
                sigma2[i]=old_sigma2;// restore old value of allele frequency
            }
            else
            {
                log_likelihood=log_likelihood+diff_log_likelihood;
                locus_likelihood[i]=new_log_likelihood_part;
                acc_sigma2[i]++;
            }
        }
    } // end of loop over loci

}
