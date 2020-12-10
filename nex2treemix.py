#!/usr/bin/env python

"""
Find biallelic SNPs in nexus files and output file of
allelic population counts for TreeMix. Will select a
single biallelic snp per locus, so snps are unlinked.
Can specify to use a random snp per locus or select snp
that has highest coverage.

date: 1 Aug 2019
author: J. Satler
version: 1

usage:
    python nex2treemix.py traits.file /path/to/nexus/files random|coverage
"""

import os
import sys
import glob
import gzip
import random
from Bio import SeqIO

def traits(infile):
    """read traits file and assign individuals to OTUs"""
    with open(infile, 'r') as t:
        return {i.split()[0]:i.split()[1] for i in t}

def locus_files(files):
    """get list of nexus files"""
    return glob.glob(files + "/*.nex*")

def read_data(data):
    """read nexus files and filter if meets sampling requirement"""
    with open(data, 'r') as locus:
        return [[ind.id, str(ind.seq)] for ind in SeqIO.parse(locus, "nexus")]

def find_biallelic_snps(mat, taxa):
    """locate snps in locus"""
    s = {n[0]:[] for n in mat}
    z_mat = [f[1].upper() for f in mat]
    coverage = []

    #find snps with allowed base pairs
    allowed = ['A', 'C', 'G', 'T']
    for index, i in enumerate(zip(*z_mat)):
        if len(set(i).intersection(allowed)) == 2:
            #found a biallelic snp
            #skip if not all pops sampled at least once
            p = {sp:0 for sp in set(taxa.values())}
            for a in range(len(mat)):
                if z_mat[a][index].upper() in allowed:
                    p[taxa[mat[a][0]]] += 1
            if all(val >= 1 for val in p.values()) == False:
                continue

            for j in range(len(mat)):
                s[mat[j][0]].append(z_mat[j][index])
            coverage.append(get_coverage(i, allowed))
    return s, coverage

def get_coverage(snp, allowed):
    """count number of allowed alleles in snp"""
    count = 0
    for allele in snp:
        if allele in allowed:
            count += 1
    return count

def get_snp(loc_snp, coverage, snp):
    """select one biallelic snp"""
    if snp == "random":
        unlinked_random_snp = random.randint(0, len(loc_snp.values()[0]) - 1)
        return {k:v[unlinked_random_snp] for (k,v) in loc_snp.items()}
    else:
        #select snp with maximum coverage
        max_coverage = max(coverage)
        position = [i for i, j in enumerate(coverage) if j == max_coverage]
        unlinked_coverage_snp = random.choice(position)
        return {k:v[unlinked_coverage_snp] for (k,v) in loc_snp.items()}

def get_pop_counts(traits, pops, snps):
    """get biallelic pop counts"""
    allowed = ['A', 'C', 'G', 'T']
    pcounts = []
    for k, v in snps.items():
        a = list(set(a for a in v.values() if a in allowed))

        #keep track of allele counts within pops
        locus_counts = {p:[0,0] for p in pops}
        for ind, allele in v.items():
            if allele == a[0]:
                locus_counts[traits[ind]][0] += 1
            elif allele == a[1]:
                locus_counts[traits[ind]][1] += 1
        pcounts.append(locus_counts)
    write_to_file(pops, pcounts)

def write_to_file(pops, pcounts):
    """write pop counts to treemix input file"""
    with gzip.open("treemixin.gz", "w") as out:
        out.write(" ".join(pops) + "\n")

        #add pop allele counts by snp
        for snp in pcounts:
            for i in pops:
                out.write(','.join(str(p) for p in snp[i]) + " ")
            out.write("\n")

def Main():
    if len(sys.argv) != 4:
        print "python nex2treemix.py traits.file /path/to/nexus/files random|coverage"
        sys.exit()

    tr = traits(sys.argv[1])
    loci = locus_files(sys.argv[2])

    pops = list(set([p for p in tr.values()]))
    #unlinked biallelic snps
    snp_loci = {}

    for index, i in enumerate(loci):
        loc = read_data(i)
        loc_snps, coverage = find_biallelic_snps(loc, tr)

        #check if locus has any snps
        if not loc_snps.values()[0]:
            continue

        if sys.argv[3].lower().startswith("r"):
            what_snp = "random"
            loc_snp_random = get_snp(loc_snps, coverage, what_snp)
            #add to master dictionary
            snp_loci[index] = loc_snp_random
        else:
            what_snp = "coverage"
            loc_snp_coverage = get_snp(loc_snps, coverage, what_snp)
            #add to master dictionary
            snp_loci[index] = loc_snp_coverage

    #get pop allele counts and write to file
    get_pop_counts(tr, pops, snp_loci)

if __name__ == '__main__':
    Main()
