# BayeScanHierachical
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

## Prior on models
Every coefficient that enters the model is penalised by the prior odds for the neutral model (`-pr_odds`, default 10). For a selection pressure that contains several groups, the models are: no selection, selection in any subset of the groups except all of them (each selected group has its own coefficient), and convergent selection, where all groups share a single coefficient. Relative to no selection, a model with k selected groups has a prior weight of `pr_odds^-k`.

The convergent model adds one coefficient and has its own prior odds, set with `-pr_odds_conv`. The default is the value of `-pr_odds`: convergent selection then has the prior weight of selection in a single group, whatever the number of groups in the pressure. With dominant markers a pressure has a single coefficient shared by its groups, and `-pr_odds_conv` applies to it when the pressure contains several groups.

Before this option existed, the convergent model of a pressure with G groups had a prior weight of `pr_odds^-(G-1)`. Nothing changes for pressures with two groups. For more groups the former behaviour is obtained with `-pr_odds_conv` set to `pr_odds^(G-1)`, for example `-pr_odds_conv 100` for three groups with the default prior odds.

The chain only reaches the convergent model from a model with all groups but one selected. With the default prior these models become much less likely than the convergent model as the number of groups grows, so the chain enters and leaves the convergent model rarely and its posterior probability is estimated with more noise. With default run lengths on the example SNP data, the typical error on this probability was 0.008 with three groups and 0.02 with four, with errors of up to 0.1 for single loci. With four groups in a pressure, use longer runs (`-n`) and compare independent runs. From five groups on, the probability of convergent selection should not be trusted without much longer runs.

A statically linked binary can be downloaded in [Releases](https://github.com/mfoll/BayeScanHierachical/releases), it should work with most Linux distributions. A docker container is available [here](https://hub.docker.com/r/mfoll/bayescanhierachical/builds/).

## Reference: 
Foll M, Gaggiotti OE, Daub JT, Vatsiou A, Excoffier L. Widespread signals of convergent adaptation to high altitude in Asia and america. Am J Hum Genet. 2014 Oct 2;95(4):394-407. doi: 10.1016/j.ajhg.2014.09.002. Epub 2014 Sep 25. PubMed PMID: [25262650](https://www.ncbi.nlm.nih.gov/pubmed/25262650); PubMed Central PMCID: [PMC4185124](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185124/). 
