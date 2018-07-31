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

## Reference: 
Foll M, Gaggiotti OE, Daub JT, Vatsiou A, Excoffier L. Widespread signals of convergent adaptation to high altitude in Asia and america. Am J Hum Genet. 2014 Oct 2;95(4):394-407. doi: 10.1016/j.ajhg.2014.09.002. Epub 2014 Sep 25. PubMed PMID: [25262650](https://www.ncbi.nlm.nih.gov/pubmed/25262650); PubMed Central PMCID: [PMC4185124](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185124/). 
