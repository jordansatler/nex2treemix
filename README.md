# nex2treemix

This script will generate an input file of allele population counts for *TreeMix* (Pickrell & Pritchard 2012). Input data are a set of nexus files and a traits file linking individual to population. The script will select a single biallelic, unlinked SNP per locus. If multiple SNPs are present within a locus, you can specify to select a SNP at random or select the SNP with the highest sampling coverage. If multiple SNPs are tied for highest sampling coverage, one of those SNPs will be selected at random.

usage:  
```python
    python nex2treemix.py traits.file /path/to/nexus/files random|coverage
```
You need to specify a traits file, path to folder of nexus files, and if you want unlinked SNPs sampled at random or sampled based on highest coverage.

***
The traits file is tab-delimited and assigns individuals to populations:

ind1a&nbsp;&nbsp;&nbsp;&nbsp;pop1  
ind1b&nbsp;&nbsp;&nbsp;&nbsp;pop1  
ind2a&nbsp;&nbsp;&nbsp;&nbsp;pop1  
ind2b&nbsp;&nbsp;&nbsp;&nbsp;pop1  
ind3a&nbsp;&nbsp;&nbsp;&nbsp;pop2  
ind3b&nbsp;&nbsp;&nbsp;&nbsp;pop2  
ind4a&nbsp;&nbsp;&nbsp;&nbsp;pop2  
ind4b&nbsp;&nbsp;&nbsp;&nbsp;pop2  
...
***
References:

Pickrell JK, Pritchard JK (2012) Inference of Population Splits and Mixtures from Genome-Wide Allele Frequency Data. PLoS Genet 8(11): e1002967. [https://doi.org/10.1371/journal.pgen.1002967](https://doi.org/10.1371/journal.pgen.1002967)
