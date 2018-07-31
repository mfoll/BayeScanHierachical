/*
This header file contains the definitions of data structures.
1) A structure that stores the characteristics of each locus, i.e.
number of alleles and allele frequency spectrum.
2) A structure that stores the characteristics of each population, i.e.
sample size and allele frequency spectrum for each locus.
*/

#include <vector>

struct indiv_data  // selection2
{
    float intensity;
};

struct locus_data
{
    int n; // total number of alleles in the sample = n_ij
    int nA1;  /* observed allele frequency spectrum = nA1_ij*/
    float p; // estimated allele frequency = p_ij
    double mean_p; // mean value of p so far (after pilots and burn-in)
    struct indiv_data *indiv; // selection2

    int ar; /*number of allelic classes = K_i*/
    int alleleCount; // total number of alleles in the sample = n_ij
    int *data_allele_count;//[max_nr_alleles];  /* observed allele frequency spectrum = a_ij*/
};

struct pop_data
{
    locus_data *locus; // alleles count and freq
    int group; // group index (between 1 and G)
};

struct group_locus_data
{
    double p; // estimated allele frequency in group, array of size I (nb loci)
    double mean_p; // mean value of p so far (after pilots and burn-in)
    int ar; /*number of allelic classes = K_i*/
    double *allele;//[max_nr_alleles]; // allele frequency
};

struct group_data
{
    group_locus_data *locus;
    std::vector<int> member; // vector of populations index member of the group
    int pressure; // pressure index between 1 and P
};

struct pressure_data
{
    std::vector<int> member; // vector of group index member sharing that pressure
};

struct allele_freq  // for ancestral population
{
    int ar; /*number of allelic classes = K_i*/
    double *allele;//[max_nr_alleles]; // allele frequency
};

