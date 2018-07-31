#include "global_defs.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "errors.cpp"
#include <algorithm>

#include <vector>


using namespace std;

int read_structure(string struct_file_name, std::ofstream& outfile) {

    ifstream infile(struct_file_name.c_str());
    assure(infile);

    string line;

    while (getline(infile, line, '=')) {
        if (line.length() >= 8 && line.substr(line.length() - 8, line.length()) == "[groups]") //read nb of groups
        {
            getline(infile, line);
            istringstream read_line(line);
            read_line >> G;
            outfile << "There are " << G << " groups." << endl << endl;

            try {
                group = new group_data[G];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for group_data " << endl;
                return 1;
            }
            for (int g = 0; g < G; g++) {
                try {
                    group[g].locus = new group_locus_data[I];
                } catch (const std::exception & Exp) {
                    cout << "Not enough memory for group[" << (g + 1) << "].locus " << endl;
                    return 1;
                }
            }

            for (int g = 0; g < G; g++) // locus index
            {
                getline(infile, line); // read the line of group g members
                istringstream read_line(line);

                int group_index; // check the group index is correct
                read_line >> group_index;
                int nb_pops;
                read_line >> nb_pops;
                for (int j = 0; j < nb_pops; j++) {
                    int cur_pop;
                    read_line >> cur_pop; // read population in group
                    group[group_index - 1].member.push_back(cur_pop - 1);
                    pop[cur_pop - 1].group = group_index - 1;
                }
            }

            // put allele number to zero, used after to initialize only once vector of allele freq
            for (int g = 0; g < G; g++) {
                for (int i = 0; i < I; i++)
                    group[g].locus[i].ar = 0;
            }
            try
            {
                eta2=new double*[G];
                for (int g=0; g < G; g++)
                    eta2[g]= new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for eta2" << endl;
                return 1;
            }
            try {
                eta2_included = new bool*[G];
                for (int g = 0; g < G; g++)
                    eta2_included[g] = new bool[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for eta2_included" << endl;
                return 1;
            }
            try {
                mean_eta2 = new double*[G];
                for (int g = 0; g < G; g++)
                    mean_eta2[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for mean_eta2" << endl;
                return 1;
            }
            try {
                var_eta2 = new double*[G];
                for (int g = 0; g < G; g++)
                    var_eta2[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for var_eta2" << endl;
                return 1;
            }
            try {
                post_eta2 = new double*[G];
                for (int g = 0; g < G; g++)
                    post_eta2[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for post_eta2" << endl;
                return 1;
            }
            try {
                nb_eta2 = new int*[G];
                for (int g = 0; g < G; g++)
                    nb_eta2[g] = new int[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for nb_eta2" << endl;
                return 1;
            }
            try {
                acc_eta2 = new double*[G];
                for (int g = 0; g < G; g++)
                    acc_eta2[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for acc_eta2" << endl;
                return 1;
            }
            try {
                var_prop_eta2 = new double*[G];
                for (int g = 0; g < G; g++)
                    var_prop_eta2[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for var_prop_eta2" << endl;
                return 1;
            }

            
            try {
                beta = new double[G];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for beta" << endl;
                return 1;
            }

            try // selection2
            {
                acc_beta = new double[G];
                var_prop_beta = new double[G];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for acc_rte beta" << endl;
                return 1;
            }
            try {
                cur_fsc = new double*[G];
                for (int g = 0; g < G; g++)
                    cur_fsc[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for cur_fsc" << endl;
                return 1;
            }
            try {
                post_fsc = new double*[G];
                for (int g = 0; g < G; g++)
                    post_fsc[g] = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for post_fsc" << endl;
                return 1;
            }
            // write out file check
            for (int g = 0; g < G; g++) {
                outfile << "There are " << group[g].member.size() << " pops in group " << (g + 1) << ": ";
                for (int k = 0; k < group[g].member.size(); k++) {
                    outfile << (group[g].member[k] + 1) << " ";
                }
                outfile << endl;
            }
            outfile << endl;

        }
        else if (line.length() >= 11 && line.substr(line.length() - 11, line.length()) == "[pressures]") //read nb of groups
        {
            getline(infile, line);
            istringstream read_line(line);
            read_line >> P;
            outfile << "There are " << P << " pressures." << endl << endl;

            try {
                pressure = new pressure_data[P];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for pressure_data " << endl;
                return 1;
            }

            // read pressure info
            for (int p = 0; p < P; p++) // locus index
            {
                getline(infile, line); // read the line of group g members
                istringstream read_line(line);

                int pressure_index; // check the group index is correct
                read_line >> pressure_index;
                int nb_groups;
                read_line >> nb_groups;
                for (int g = 0; g < nb_groups; g++) {
                    int cur_group;
                    read_line >> cur_group; // read population in group
                    pressure[pressure_index - 1].member.push_back(cur_group - 1);
                    group[cur_group - 1].pressure = pressure_index - 1;
                }
            }

            // write check file
            for (int p = 0; p < P; p++) {
                outfile << "There are " << pressure[p].member.size() << " groups in pressure " << (p + 1) << ": ";
                for (int k = 0; k < pressure[p].member.size(); k++) {
                    outfile << (pressure[p].member[k] + 1) << ", ";
                }
                outfile << endl;
            }
            outfile << endl;
        }
    }

    // if no pressure info, I automatically create one pressure per group
    if (P == 0) {
        P = G;
        try {
            pressure = new pressure_data[P];
        } catch (const std::exception & Exp) {
            cout << "Not enough memory for pressure_data " << endl;
            return 1;
        }

        // read pressure info
        for (int p = 0; p < P; p++) // locus index
        {
            pressure[p].member.push_back(p);
            group[p].pressure = p;
        }

    }

    try {
        eta = new double*[P];
        for (int p = 0; p < P; p++)
            eta[p] = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for eta" << endl;
        return 1;
    }
    try {
        eta_included = new bool*[P];
        for (int p = 0; p < P; p++)
            eta_included[p] = new bool[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for eta_included" << endl;
        return 1;
    }
    try {
        mean_eta = new double*[P];
        for (int p = 0; p < P; p++)
            mean_eta[p] = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for mean_eta" << endl;
        return 1;
    }
    try {
        var_eta = new double*[P];
        for (int p = 0; p < P; p++)
            var_eta[p] = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for var_eta" << endl;
        return 1;
    }
    try {
        post_eta = new double*[P];
        for (int p = 0; p < P; p++)
            post_eta[p] = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for post_eta" << endl;
        return 1;
    }
    try {
        nb_eta = new int*[P];
        for (int p = 0; p < P; p++)
            nb_eta[p] = new int[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for nb_eta" << endl;
        return 1;
    }
    try {
        acc_eta = new double*[P];
        for (int p = 0; p < P; p++)
            acc_eta[p] = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_eta" << endl;
        return 1;
    }
    try {
        var_prop_eta = new double*[P];
        for (int p = 0; p < P; p++)
            var_prop_eta[p] = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for var_prop_eta" << endl;
        return 1;
    }

    infile.close();
    infile.clear();

}

int read_input(std::ifstream& infile, string struct_file_name, string outcheck) {

    ofstream outfile;


    outfile.open(outcheck.c_str());
    assure(outfile);

    outfile << "Summary of parameters and input files." << endl;
    outfile << "Please check that all is correct while calculation is starting..." << endl << endl;


    string line;

    while (getline(infile, line, '=')) {
        if (line.length() >= 6 && line.substr(line.length() - 6, line.length()) == "[loci]") //read nb of loci
        {
            getline(infile, line);
            istringstream read_line(line);
            read_line >> I;
            outfile << "There are " << I << " loci." << endl << endl;
            try {
                alpha = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for alpha" << endl;
                return 1;
            }
            try {
                alpha_included = new bool[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for alpha_included" << endl;
                return 1;
            }
            try {
                discarded_loci = new bool[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for discarded_loci" << endl;
                return 1;
            }
            try {
                mean_alpha = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for mean_alpha" << endl;
                return 1;
            }
            try {
                var_alpha = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for var_alpha" << endl;
                return 1;
            }
            try {
                post_alpha = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for post_alpha" << endl;
                return 1;
            }
            /*try {
                  p_value=new double[I];
            }
            catch ( const std::exception & Exp ) {
                cout << "Not enough memory for p_value"<< endl;
                return 1;
            } */
            try {
                nb_alpha = new int[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for nb_alpha" << endl;
                return 1;
            }
            try {
                cur_fct = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for cur_fct" << endl;
                return 1;
            }
            try {
                post_fct = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for post_fct" << endl;
                return 1;
            }
            try {
                freq_ancestral = new double[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for freq_ancestral" << endl;
                return 1;
            }
            if (codominant == 1) // selection2 : I forgot this before (waste of unused memory) ?
            {
                try {
                    freq_locus = new allele_freq[I];
                } catch (const std::exception & Exp) {
                    cout << "Not enough memory for allele_freq " << endl;
                    return 1;
                }
            }

            //selection2
            try {
                acc_alpha = new double [I];
                var_prop_alpha = new double [I];
                acc_freq_ancestral = new double [I];
                e_ancestral = new double [I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for acc_rte alpha or freq_ancestral " << endl;
                return 1;
            }

        } else if (line.length() >= 13 && line.substr(line.length() - 13, line.length()) == "[populations]") //read nb of populations
        {
            getline(infile, line);
            istringstream read_line(line);
            read_line >> J;
            outfile << "There are " << J << " populations." << endl << endl;
            try {
                pop = new pop_data[J];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for pop_data " << endl;
                return 1;
            }

            try {
                theta = new double[J];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for theta" << endl;
                return 1;
            }

            try // selection2
            {
                acc_theta = new double[J];
                var_prop_theta = new double[J];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for acc_rte theta" << endl;
                return 1;
            }

            if (codominant < 1) {
                try {
                    f = new double[J];
                } catch (const std::exception & Exp) {
                    cout << "Not enough memory for f" << endl;
                    return 1;
                }
            }

            if (codominant < 1) // selection2
            {
                try {
                    acc_f = new double[J];
                    e_f = new double[J];
                } catch (const std::exception & Exp) {
                    cout << "Not enough memory for acc_rte f" << endl;
                    return 1;
                }
            }

        }

        else if (line.length() >= 5 && line.substr(line.length() - 5, line.length()) == "[pop]") //read allele count
        {
            // read population structure if first time
            if (P==0)
                read_structure(struct_file_name, outfile);

            // initialize fst table
            /*fst=(double**)malloc(I*sizeof(double *));
            for(int i=0;i<I;i++)
            fst[i]=(double*)malloc(J*sizeof(double));*/

            try {
                acc_freq = new double*[I];
                e_freq = new double*[I];
                for (int i = 0; i < I; i++) {
                    acc_freq[i] = new double[J];
                    e_freq[i] = new double[J];
                }
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for acc_freq" << endl;
                return 1;
            }

            try {
                acc_freq_group = new double*[I];
                e_freq_group = new double*[I];
                for (int i = 0; i < I; i++) {
                    acc_freq_group[i] = new double[G];
                    e_freq_group[i] = new double[G];
                }
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for acc_freq_group" << endl;
                return 1;
            }

            getline(infile, line);
            istringstream read_line(line);
            int j;
            read_line >> j; // population index
            j = j - 1;

            try {
                pop[j].locus = new locus_data[I];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for pop[" << (j + 1) << "].locus " << endl;
                return 1;
            }

            for (int i = 0; i < I; i++) // locus index
            {
                getline(infile, line); // read the line of allele count at locus i
                istringstream read_line(line);

                int locus_index; // check the locus index is correct
                read_line >> locus_index;
                if (locus_index != i + 1) {
                    cout << "Could not find allele count for population " << (j + 1) << " at locus " << (i + 1) << endl;
                    return 1;
                }

                if (codominant == 1) // selection2
                {
                    read_line >> pop[j].locus[i].alleleCount; // total number of alleles in the sample at locus i

                    read_line >> pop[j].locus[i].ar; // number of allelic classes at locus i
                    
                    bool ar1=false;
                    if (pop[j].locus[i].ar==1) {
						ar1=true;
						pop[j].locus[i].ar=2;
                    }
                    
                    freq_locus[i].ar = pop[j].locus[i].ar; // the number is the same for all populations

                    try {
                        pop[j].locus[i].data_allele_count = new int[pop[j].locus[i].ar];
                    } catch (const std::exception & Exp) {
                        cout << "Not enough memory for pop[" << (j + 1) << "locus[" << (i + 1) << "].data_allele_count" << endl;
                        return 1;
                    }
                    if (j == 0) {
                        try {
                            freq_locus[i].allele = new double[freq_locus[i].ar];
                        } catch (const std::exception & Exp) {
                            cout << "Not enough memory for freq_locus[" << (i + 1) << "].allele" << endl;
                            return 1;
                        }
                    }
                    if (group[pop[j].group].locus[i].ar == 0) // initialise for groups only if not done already
                    {
                        group[pop[j].group].locus[i].ar = pop[j].locus[i].ar; // and also for groups
                        try {
                            group[pop[j].group].locus[i].allele = new double[group[pop[j].group].locus[i].ar];
                        } catch (const std::exception & Exp) {
                            cout << "Not enough memory for group[" << (pop[j].group + 1) << "].locus[" << (i + 1) << "].allele" << endl;
                            return 1;
                        }
                    }

					if (!ar1)
                    {
                    	for (int k = 0; k < pop[j].locus[i].ar; k++)
                        	read_line >> pop[j].locus[i].data_allele_count[k];
                    }
                    else
                    {
						read_line >> pop[j].locus[i].data_allele_count[0];
						pop[j].locus[i].data_allele_count[1]=0;
		    		}

                } else {
                    read_line >> pop[j].locus[i].n;
                    read_line >> pop[j].locus[i].nA1;
                }

            }

        }

    }

    outfile << "Burn in: " << discard << endl;

    outfile << "Thining interval: " << interval << endl;

    outfile << "Sample size: " << nr_out << endl;

    outfile << "Resulting total number of iterations: " << tot_nr_of_iter << endl;

    outfile << "Nb of pilot runs: " << nb_pilot << endl;

    outfile << "Length of each pilot run: " << pilot_length << endl << endl;


    // write allele counts
    outfile << "Allele counts:" << endl;
    for (int j = 0; j < J; j++) {
        for (int i = 0; i < I; i++) {
            outfile << "Pop. " << j + 1 << " locus " << i + 1 << " : ";
            if (codominant == 1) // selection2
            {
                for (int k = 0; k < pop[j].locus[i].ar; k++)
                    outfile << setw(3) << pop[j].locus[i].data_allele_count[k] << " ";
            } else
                outfile << setw(3) << pop[j].locus[i].n << " " << pop[j].locus[i].nA1;
            outfile << endl;
        }
        outfile << endl;
    }
    outfile.close();

    infile.close();

    return 0;
}

int read_input_intensity(string outname, string struct_file_name, string outcheck) {
    ofstream outfile;

    outfile.open(outcheck.c_str());
    assure(outfile);

    outfile << "Summary of parameters and input files." << endl;
    outfile << "Please check that all is correct while calculation is starting..." << endl << endl;

    vector<int> individuals; // size should be always the nb of pop , contains the number of indiv in each pop
    double trash;
    int cur_pop;
    int cur_indiv;
    //    int cur_group;

    // first pass to read number of pop, loci, and indivs
    ifstream infile(outname.c_str());
    assure(infile);

    I = 0;
    J = 0;
    //    G = 0;

    string line; // contains each line

    // use the first line for number of loci
    getline(infile, line);
    istringstream read_line(line);
    read_line >> cur_indiv;
    read_line >> cur_pop;
    //    read_line >> cur_group;
    if (cur_pop > J) {
        J = cur_pop;
        individuals.resize(J, 0);
    }
    //    if (cur_group > G)
    //        G = cur_group;
    if (cur_indiv > individuals[cur_pop - 1])
        individuals[cur_pop - 1] = cur_indiv;

    //while(read_line.good())
    while (read_line >> trash) {
        //	read_line >> trash;
        I++;
    }

    //I--;

    while (getline(infile, line)) {
        istringstream read_line(line);
        read_line >> cur_indiv;
        read_line >> cur_pop;
        //        read_line >> cur_group;
        if (cur_pop > J) {
            J = cur_pop;
            individuals.resize(J, 0);
        }
        //        if (cur_group > G)
        //            G = cur_group;
        if (cur_indiv > individuals[cur_pop - 1])
            individuals[cur_pop - 1] = cur_indiv;
    }


    infile.close();
    infile.clear();

    outfile << "There are " << I << " loci." << endl << endl;

    outfile << "There are " << J << " populations." << endl << endl;

    // read structure file
    read_structure(struct_file_name, outfile);

    // second pass to fill the database
    infile.open(outname.c_str());

    // allocate memory
    try {
        locus_likelihood = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for locus_likelihood" << endl;
        return 1;
    }
    try {
        alpha = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for alpha" << endl;
        return 1;
    }
    try {
        discarded_loci = new bool[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for discarded_loci" << endl;
        return 1;
    }
    try {
        alpha_included = new bool[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for alpha_included" << endl;
        return 1;
    }
    try {
        mean_alpha = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for mean_alpha" << endl;
        return 1;
    }
    try {
        var_alpha = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for var_alpha" << endl;
        return 1;
    }
    try {
        post_alpha = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for post_alpha" << endl;
        return 1;
    }
    try {
        nb_alpha = new int[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for nb_alpha" << endl;
        return 1;
    }
    try {
        cur_fct = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for cur_fct" << endl;
        return 1;
    }
    try {
        post_fct = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for post_fct" << endl;
        return 1;
    }
    try {
        freq_ancestral = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for freq_ancestral" << endl;
        return 1;
    }
    try {
        mu = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for mu" << endl;
        return 1;
    }
    try {
        delta = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for delta" << endl;
        return 1;
    }
    try {
        sigma1 = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for sigma1" << endl;
        return 1;
    }
    try {
        sigma2 = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for sigma2" << endl;
        return 1;
    }
    try {
        max_intenstiy = new double[I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for max_intenstiy" << endl;
        return 1;
    }

    try {
        acc_alpha = new double [I];
        var_prop_alpha = new double [I];
        acc_freq_ancestral = new double [I];
        e_ancestral = new double [I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_rte alpha or freq_ancestral " << endl;
        return 1;
    }

    try {
        acc_mu = new double [I];
        acc_delta = new double [I];
        acc_sigma1 = new double [I];
        acc_sigma2 = new double [I];
        e_mu = new double [I];
        e_delta = new double [I];
        var_prop_sigma1 = new double [I];
        var_prop_sigma2 = new double [I];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_rte mu, delta or sigma " << endl;
        return 1;
    }

    try {
        pop = new pop_data[J];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for pop_data " << endl;
        return 1;
    }
    try {
        f = new double[J];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for f" << endl;
        return 1;
    }

    try {
        theta = new double[J];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for theta" << endl;
        return 1;
    }

    try // selection2
    {
        acc_theta = new double[J];
        var_prop_theta = new double[J];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_rte theta" << endl;
        return 1;
    }


    try {
        acc_f = new double[J];
        e_f = new double[J];
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_rte f" << endl;
        return 1;
    }

    for (int j = 0; j < J; j++) {
        try {
            pop[j].locus = new locus_data[I];
        } catch (const std::exception & Exp) {
            cout << "Not enough memory for pop[" << (j + 1) << "].locus ";
            return 1;
        }
    }

    try {
        acc_freq = new double*[I];
        e_freq = new double*[I];
        for (int i = 0; i < I; i++) {
            acc_freq[i] = new double[J];
            e_freq[i] = new double[J];
        }
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_freq" << endl;
        return 1;
    }

    for (int j = 0; j < J; j++) {
        for (int i = 0; i < I; i++) {
            try {
                pop[j].locus[i].indiv = new indiv_data[individuals[j]];
            } catch (const std::exception & Exp) {
                cout << "Not enough memory for pop[" << (j + 1) << "].locus ";
                return 1;
            }
            pop[j].locus[i].n = individuals[j];
        }
    }

    for (int i = 0; i < I; i++)
        max_intenstiy[i] = 0;


    try {
        acc_freq_group = new double*[I];
        e_freq_group = new double*[I];
        for (int i = 0; i < I; i++) {
            acc_freq_group[i] = new double[G];
            e_freq_group[i] = new double[G];
        }
    } catch (const std::exception & Exp) {
        cout << "Not enough memory for acc_freq_group" << endl;
        return 1;
    }

    while (getline(infile, line)) {
        istringstream read_line(line);
        read_line >> cur_indiv;
        read_line >> cur_pop;
        //        read_line >> cur_group;

        //        group[cur_group - 1].member.push_back(cur_pop - 1);
        //        pop[cur_pop - 1].group = cur_group - 1;

        // read intensities
        for (int i = 0; i < I; i++) {
            read_line >> pop[cur_pop - 1].locus[i].indiv[cur_indiv - 1].intensity;
            if (pop[cur_pop - 1].locus[i].indiv[cur_indiv - 1].intensity > max_intenstiy[i])
                max_intenstiy[i] = pop[cur_pop - 1].locus[i].indiv[cur_indiv - 1].intensity;
        }
    }

    if (!SNP_genotypes) {
        for (int i = 0; i < I; i++) {
            for (int j = 0; j < J; j++) {
                for (int k = 0; k < pop[j].locus[i].n; k++) {
                    if (max_intenstiy[i] > 0)
                        pop[j].locus[i].indiv[k].intensity /= max_intenstiy[i];
                }
            }
        }
    }

    // make vectors of pop membership in groups unique
//    for (int g = 0; g < G; g++) {
//        std::sort(group[g].member.begin(), group[g].member.end());
//        group[g].member.erase(unique(group[g].member.begin(), group[g].member.end()), group[g].member.end());
//    }

    infile.close();
    infile.clear();

    // write check file
    for (int j = 0; j < J; j++) {
        outfile << "There are " << individuals[j] << " individuals in population " << j + 1 << "." << endl;
    }
    outfile << endl;

    outfile << "Burn in: " << discard << endl;

    outfile << "Thining interval: " << interval << endl;

    outfile << "Sample size: " << nr_out << endl;

    outfile << "Resulting total number of iterations: " << tot_nr_of_iter << endl;

    outfile << "Nb of pilot runs: " << nb_pilot << endl;

    outfile << "Length of each pilot run: " << pilot_length << endl << endl;

    outfile.close();


    return 0;
}

//////////////////////////////
// read discarded loci file
//////////////////////////////

int read_discarded(std::ifstream& infile) {
    int locus;
    while (infile >> locus) {
        if (locus >= 1 && locus <= I)
            discarded_loci[locus - 1] = true;
    }
}

///////////////////////////////////
// write output file
///////////////////////////////////

void write_output(std::ofstream& outfile) {
    // write iteration index
    outfile << iter << "  ";
    // write loglikelihood
    outfile << setprecision(8) << " " << log_likelihood;
    if (codominant < 1) {
        // write a
        //outfile << setprecision(8) << " " << a_p;
        // write fis
        for (int j = 0; j < J; j++)
            outfile << setprecision(8) << " " << f[j];
    }
    // write fct
    for (int g = 0; g < G; g++)
        outfile << setprecision(8) << " " << 1 / (1 + exp(-beta[g]));
    // write fsc
    for (int j = 0; j < J; j++)
        outfile << setprecision(8) << " " << 1 / (1 + exp(-theta[j]));
    // write prop of non zero alpha
    //outfile << setprecision(8) << " " << ((double)nb_alpha_included)/((double)I);
    // write alpha
    if (!fstat  && all_trace) {
        for (int i = 0; i < I; i++) {
            if (!discarded_loci[i])
                outfile << setprecision(8) << " " << alpha[i];
        }
        for (int p = 0; p < P; p++) {
            for (int i = 0; i < I; i++) {
                if (!discarded_loci[i])
                    outfile << setprecision(8) << " " << eta[p][i];
            }
        }
        for (int g = 0; g < G; g++) {
            for (int i = 0; i < I; i++) {
                if (!discarded_loci[i])
                    outfile << setprecision(8) << " " << eta2[g][i];
            }
        }
    }
    // write alpha_included
    //for (int i=0;i<I;i++)
    //outfile << setprecision(8) << " " << (int)alpha_included[i];
    // write fst
    //for (int i=0;i<I;i++)
    //outfile << setprecision(8) << " " << cur_fct[i];

    outfile << endl;
}


////////////////////////////////////////////////////////////////
// write output file containing models names in decimal format
//////////////////////////////////////////////////////////////
// ...dcba: a is for alpha. b for eta_0, c for eta_1 etc.
// ...0000: neutral model
// ...0001: alpha included
// ...0101: eta_1 and alpha included etc.

void write_models(std::ofstream& outfile_models) {
    outfile_models << iter << "  ";
    for (int i = 0; i < I; i++) {
        if (!discarded_loci[i]) {
            int m = (alpha_included[i]);
            for (int p = 0; p < P; p++)
                m += (int) pow(2, p + 1)*(eta_included[p][i]);
            outfile_models << m << " ";
        }
    }
    outfile_models << endl;
}

///////////////////////////////////
// write output file for mixture model of AFLP band intensity
///////////////////////////////////

void write_output_intensity(std::ofstream& outfile) // selection2
{
    // write iteration index
    outfile << iter << "  ";
    // write mu and mu+delta
    for (int i = 0; i < I; i++) {
        if (!discarded_loci[i]) {
            if (max_intenstiy[i] != 0)
                outfile << setprecision(8) << " " << mu[i] * max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << mu[i];

            if (max_intenstiy[i] != 0)
                outfile << setprecision(8) << " " << (delta[i] + mu[i]) * max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << delta[i] + mu[i];
        }
    }
    // write sigma1 and sigma2
    for (int i = 0; i < I; i++) {
        if (!discarded_loci[i]) {
            if (max_intenstiy[i] != 0)
                outfile << setprecision(8) << " " << sigma1[i] * max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << sigma1[i];
            if (max_intenstiy[i] != 0)
                outfile << setprecision(8) << " " << sigma2[i] * max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << sigma2[i];
        }
    }
    outfile << endl;
}


///////////////////////////////////
// write genotype freq Aa / AA for mixture model of AFLP band intensity
///////////////////////////////////

void write_output_genotype(std::ofstream outfile[]) // selection2
{
    double p, pAa, pAA;
    for (int j = 0; j < J; j++) {
        // write iteration index
        outfile[j] << iter << "  ";
        // write pAa and pAA
        for (int i = 0; i < I; i++) {
            if (!discarded_loci[i]) {
                p = pop[j].locus[i].p;
                pAA = p * p + f[j] * p * (1 - p);
                pAa = 2 * p * (1 - p)*(1 - f[j]);
                outfile[j] << setprecision(8) << " " << pAa << " " << pAA;
            }
        }
        outfile[j] << endl;
    }
}

/////////////////////////////////////////////
// write output file for allele frequencies
/////////////////////////////////////////////

void write_freq(std::ofstream outfiles[]) {
    for (int j = 0; j < J; j++) {
        outfiles[j] << iter;
        for (int i = 0; i < I; i++) //cycle over loci
        {
            if (!discarded_loci[i])
                outfiles[j] << " " << pop[j].locus[i].p;
        }
        outfiles[j] << endl;
    }
}

/////////////////////////////////////////////
// write output file for allele frequencies in groups
/////////////////////////////////////////////

void write_group_freq(std::ofstream outfiles[]) {
    for (int g = 0; g < G; g++) {
        outfiles[g] << iter;
        for (int i = 0; i < I; i++) //cycle over loci
        {
            if (!discarded_loci[i]) {
                if (codominant == 1) // selection2
                {
                    outfiles[g] << " " << group[g].locus[i].ar;
                    for (int k = 0; k < group[g].locus[i].ar; k++)
                        outfiles[g] << " " << group[g].locus[i].allele[k];
                } else
                    outfiles[g] << " " << group[g].locus[i].p;
            }
        }
        outfiles[g] << endl;
    }
}



/////////////////////////////////////////////
// write single output file for posterior allele frequencies
/////////////////////////////////////////////
// public version: added "popxx" in front of each line

void write_freq_matrix(std::ofstream& freq_pop, /*std::ofstream& freq_group,*/int cur_out) {
    // write header
    for (int i = 0; i < I; i++) {
        if (!discarded_loci[i])
            freq_pop << "locus" << (i + 1) << " ";
    }
    freq_pop << endl;
    for (int j = 0; j < J; j++) {
        freq_pop << "pop" << (j + 1) << " ";
        for (int i = 0; i < I; i++) {
            if (!discarded_loci[i])
                freq_pop << pop[j].locus[i].mean_p / ((double) cur_out) << " ";
        }
        freq_pop << endl;
    }

    // write header
    /*    for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                freq_group << "locus" << (i+1) << " ";
        }
        freq_group << endl;
        for (int g=0;g<G;g++)
        {
            freq_group << (g+1) << " ";
            for (int i=0;i<I;i++)
            {
                if (!discarded_loci[i])
                    freq_group << group[g].locus[i].mean_p/((double)cur_out) << " ";
            }
            freq_group << endl;
        }*/
}

/////////////////////////////////////////////
// write output file for ancestral allele frequencies
/////////////////////////////////////////////

void write_anc_freq(std::ofstream& outfile) {
    outfile << iter;

    for (int i = 0; i < I; i++) //cycle over loci
    {
        if (!discarded_loci[i]) {
            if (codominant == 1) // selection2
            {
                outfile << " " << freq_locus[i].ar;
                for (int k = 0; k < freq_locus[i].ar; k++)
                    outfile << " " << freq_locus[i].allele[k];
            } else
                outfile << " " << freq_ancestral[i];
        }
    }
    outfile << endl;
}

