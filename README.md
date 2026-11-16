# BayeScanHierachical

## Important: errors corrected in version 1.1 (November 2026)

**Version 1.0 contains two errors that make the posterior probabilities of selection too high. Please use version 1.1 or later, and run again the analyses made with version 1.0.**

- In the reversible-jump steps, the moves that remove a selection coefficient, except the one that leaves the convergent model, computed their acceptance probability after setting the coefficient to zero, and were accepted too rarely. The posterior probabilities of selection were overestimated and the q-values underestimated, mostly for loci with intermediate evidence. This is the error corrected in [BayeScan 2.2](https://github.com/mfoll/BayeScan).
- For a selection pressure that contains three groups or more, one proposal ratio was wrong by a factor of two. The models with all groups but one selected, and the convergent model, had twice their prior weight. Pressures with two groups are not concerned.

As an example, take the SNP example file of BayeScan (`test_SNPs.txt`) with the default options. With its ten populations in five groups of two and one pressure per group, 24 loci have a q-value below 0.05 for selection in at least one group with version 1.0, and 12 with version 1.1. With two groups of five populations and one pressure shared by the two groups, 13 loci have a q-value below 0.05 for convergent selection with both versions. Running longer chains does not correct the errors.

The errors are in the code and not in the model.

Version 1.1 also changes the prior of the convergent model for pressures with three groups or more, and adds a move to the chain. Neither is a correction, and both are described in *Prior on models* below. With three groups or more in a pressure, results therefore differ from those of version 1.0 for this reason as well; the option `-pr_odds_conv` gives the former prior.

Many thanks to Luis Mijangos and his colleagues of the [dartR](https://github.com/green-striped-gecko/dartR) team, who found the first error in BayeScan and contributed its correction there ([pull request #9 of BayeScan](https://github.com/mfoll/BayeScan/pull/9)). The same error, and the second one, were then found here.

## Usage
Hierachical version of BayeScan to detect natural selection from population-based genetic data. Most of the code and options are taken from the [original BayeScan](https://github.com/mfoll/BayeScan). There is one additional input file used to describe the hierachical structure of the populations (`-s` option). One first need to indicate the groups of populations, for example:
```
[groups]=2
1 2 1 2
2 2 3 4
```
indicates that we have two groups of population. The first group (`1`) is composed of two populations (`2`) that are labeled as `1` and `2` in the main input file. The second group (`2`) is composed of two populations (`2`) that are labeled as `3` and `4` in the main input file.

The second part of this input file describes the pattern of convergent selection tested by assigning common selection pressures to the groups of populations previously defined. For example if one wants to test for converget selection in the two groups of population:
```
[pressures]=1
1 2 1 2 
```
indicates that there is one selection pressure containing two groups labeled as `1` and `2`. When convergent selection if not of a particular interest, the pressures can be defined as:
```
[pressures]=2
1 1 1 
2 1 2
```
indicating an independant selection pressure in each group of populations.

The main input file has the format of BayeScan for codominant markers or for binary dominant markers. Band-intensity AFLP data and SNP genotype matrices (`-snp`) are not supported: the program crashes when it reads them, in version 1.0 as in 1.1.

## Prior on models
Every coefficient that enters the model is penalised by the prior odds for the neutral model (`-pr_odds`, default 10). For a selection pressure that contains several groups, the models are: no selection, selection in any subset of the groups except all of them (each selected group has its own coefficient), and convergent selection, where all groups share a single coefficient. Relative to no selection, a model with k selected groups has a prior weight of `pr_odds^-k`.

The convergent model adds one coefficient and has its own prior odds, set with `-pr_odds_conv`. The default is the value of `-pr_odds`: convergent selection then has the prior weight of selection in a single group, whatever the number of groups in the pressure. With dominant markers a pressure has a single coefficient shared by its groups, and `-pr_odds_conv` applies to it when the pressure contains several groups.

Before this option existed, the convergent model of a pressure with G groups had a prior weight of `pr_odds^-(G-1)`. Nothing changes for pressures with two groups. For more groups the former behaviour is obtained with `-pr_odds_conv` set to `pr_odds^(G-1)`, for example `-pr_odds_conv 100` for three groups with the default prior odds.

The chain reaches the convergent model either directly from the model without selection or from a model with all groups but one selected. With default run lengths on the example SNP data and two to five groups in a pressure, the typical error on the posterior probability of convergent selection was 0.004 to 0.006, and at most 0.08 for a single locus with five groups.

Executables for Linux (statically linked, it should work with most distributions) and macOS (Apple Silicon) can be downloaded in [Releases](https://github.com/mfoll/BayeScanHierachical/releases). The Docker image on [Docker Hub](https://hub.docker.com/r/mfoll/bayescanhierachical/) has not been rebuilt: it contains version 1.0 and its errors.

## Reference: 
Foll M, Gaggiotti OE, Daub JT, Vatsiou A, Excoffier L. Widespread signals of convergent adaptation to high altitude in Asia and america. Am J Hum Genet. 2014 Oct 2;95(4):394-407. doi: 10.1016/j.ajhg.2014.09.002. Epub 2014 Sep 25. PubMed PMID: [25262650](https://www.ncbi.nlm.nih.gov/pubmed/25262650); PubMed Central PMCID: [PMC4185124](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185124/). 
