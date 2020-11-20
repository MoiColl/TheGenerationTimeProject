<a name="H"></a>
# Accumulation and mutation spectrum of derived alleles

In this notebook it is shown how to reproduce the results in (X). First, polymorphic sites are called in the Simons Genome Diversity Project (SGDP) dataset ([Mallick S et al, 2016](https://doi.org/10.1038/nature18964)) using [`cpoly`](https://github.com/DReichLab/cTools) program. Then, the different filters will be applied (see below) to get polymorphic sites that were accumulated likely after the Out of Africa and outside archaic introgressed regions. Finally, the results and plots from (X) are displayed. 

As example, the crhomsome 22 is generaly used throught the notebook. The same methods can be applied for the rest of the genome.

## Contents

1. [ Library and packages ](#Lib)
2. [ Polymorphic sites per chromosome in SGDP with cpoly ](#SNPs)
3. [ Masked regions ](#Mask)
4. [ Derived allele counts per mutation type ](#Der)
    - 4.1. [ Load reference genomes ](#Ref)
    - 4.2. [ ind, snp and geno files ](#Dat)
    - 4.3. [ Filters ](#Fil)
    - 4.4. [ Mutation types per individual ](#Mut)
5. [ References ](#Bib)

<a name="Lib"></a>
## 1. Library and packages

Python packages


```python
from collections import defaultdict
import numpy as np

%load_ext rpy2.ipython
```


```r
%%R

library(tidyverse)
library(ggplot2)
```

    R[write to console]: -- Attaching packages --------------------------------------- tidyverse 1.2.1 --
    
    R[write to console]: v ggplot2 3.3.2     v purrr   0.3.4
    v tibble  3.0.1     v dplyr   1.0.0
    v tidyr   1.1.0     v stringr 1.4.0
    v readr   1.3.1     v forcats 0.5.0
    
    R[write to console]: -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    x dplyr::filter() masks stats::filter()
    x dplyr::lag()    masks stats::lag()
    


[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="SNPs"></a>
## 2. Polymorphic sites per chromosome in SGDP with cpoly

The polymorphic sites of the SGDP data ([Mallick S et al, 2016](https://doi.org/10.1038/nature18964)) are pulled using [`cpoly`](https://github.com/DReichLab/cTools), which is part of a set of tools (Ctools) to handle SGDP data, released with in the same publication.

The program requires two files to run:

1. `.par` file: The parameter file that indicates specific parameters.
2. `.ind` file: The individual file that indicates which individuals the polymorphic sites are going to be pulled from.

[`SGDP_22.par`](SGDP_22.par) and [`SGDP.ind`](SGDP.ind) files are examples used to run `cpoly` on chromosome 22 for a subset of all individuals in SGDP.


```bash
%%bash

cat SGDP_22.par
```

    indivname:	SGDP.ind
    indivoutname:	SGDP_22.ind
    snpoutname:	SGDP_22.snp
    genooutname:	SGDP_22.geno
    dbhetfa:	hetfa.dblist
    dbmask:	mask.dblist
    chrom:	22
    minfilterval:	1
    polarize:	Href



```bash
%%bash

head SGDP.ind
```

    Href	M	hg19ref
    S_Abkhasian-1	M	WestEurasia
    S_Abkhasian-2	M	WestEurasia
    S_Adygei-1	M	WestEurasia
    S_Adygei-2	F	WestEurasia
    S_Albanian-1	F	WestEurasia
    S_Aleut-1	M	CentralAsiaSiberia
    S_Aleut-2	F	CentralAsiaSiberia
    S_Altaian-1	M	CentralAsiaSiberia
    S_Ami-1	M	EastAsia


`cpoly` was run with the following comand:


$ cpoly -p SGDP_22.par
The program outputs 3 files:

1. `.ind` file: contains the list of individuals used by cpoly. Number of indidivuals = $i$.
2. `.snp` file: contains the list of polymorphic sites found by cpoly and it's basic information (position, reference allele, alternative allele). Number of indidivuals = $s$.
3. `.geno` file: matrix of dimentions ($s$, $i$) with the genotypes (0 : alternative allele homozygote, 1 : heterozygote, 2: reference allele homozygote, 9 : noncallable) of each individual, with the same order as indicated in the `.ind`, for each polymorphic site, with the same order as indicated in the `.snp`. 



```bash
%%bash

head SGDP_22.ind
```

    Href	M	hg19ref
    S_Abkhasian-1	M	WestEurasia
    S_Abkhasian-2	M	WestEurasia
    S_Adygei-1	M	WestEurasia
    S_Adygei-2	F	WestEurasia
    S_Albanian-1	F	WestEurasia
    S_Aleut-1	M	CentralAsiaSiberia
    S_Aleut-2	F	CentralAsiaSiberia
    S_Altaian-1	M	CentralAsiaSiberia
    S_Ami-1	M	EastAsia



```bash
%%bash

head SGDP_22.snp
```

     X:22_16855028   22 0     16855028 G T
     X:22_16855618   22 0     16855618 G A
     X:22_16959041   22 0     16959041 G C
     X:22_16959043   22 0     16959043 G T
     X:22_16959063   22 0     16959063 G C
     X:22_16959070   22 0     16959070 G T
     X:22_16959133   22 0     16959133 C T
     X:22_17003679   22 0     17003679 A G
     X:22_17003695   22 0     17003695 C T
     X:22_17015603   22 0     17015603 G A



```bash
%%bash

head -n 3 SGDP_22.geno
```

    2222222222222222222222222222222222222222222222222222922992222222292222922222222222222222222292222222922922222222229292222222222222222222222222292222229929222222222222222222222222222222222222222222922292229222292222212222299222222922222
    2111100221102101221122100012100111111212020200212012921991212212292122022002112111110112111121212292291921212111221291222911012122211212211201111111221191022010291122212021211220092011210101011211122911291111211020022222922229101292222
    2222229292999299999222929299299222929922292922992929299929229922929922299229922222229229222292229922292299922292292992922222922229222922922929229992929992299992922292222999229929299222299222922999292992291922992912229922292299229992229


[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="Mask"></a>
## 3. Masked regions

In this section, the February 2009 human reference sequence (GRCh37) will be obtained with repeatitive regions and archaic introgressed regions masked.

First, we download a reference sequence for each chromosome from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ in which the following is noted:

*Repeats from RepeatMasker and Tandem Repeats Finder (with period of 12 or less) are shown in lower case; non-repeating sequence is shown in upper case. RepeatMasker was run with the -s (sensitive) setting. Using: Jan 29 2009 (open-3-2-7) version of RepeatMasker and RELEASE 20090120 of library RepeatMaskerLib.embl*


The following command was used to download the reference chromosome 22:
$ rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz .
$ gunzip chr22.fa.gz
We then obtained the regions with evidence of archaic introgression by merging all fragments from chromosome 22 in all individuals used in this section (check [SGDP.ind](SGDP.ind)) from `Data1_archaicfragments.txt` with the following [bedtools](https://bedtools.readthedocs.io/en/latest/) (v2.29.2) ([Aaron R et al, 2010](https://doi.org/10.1093/bioinformatics/btq033)) command:


```bash
%%bash

awk '{if(NR > 1 && $2 != "WHG" && $2 != "LBK" && $2 != "Ust_ishim" && $3 == 22){print "chr"$3"\t"$4"\t"$5}}' Data1_archaicfragments.txt | sort -k1n -k2n -k3n | bedtools merge -i - | head
```

    chr22	17257000	17443000
    chr22	17468000	17676000
    chr22	17679000	17761000
    chr22	17782000	18082000
    chr22	18250000	18263000
    chr22	18295000	18379000
    chr22	18434000	18452000
    chr22	18458000	18484000
    chr22	18499000	18523000
    chr22	18524000	18534000


and masked those regions from the reference chromosome 22 fasta file:


```bash
%%bash

bedtools maskfasta -fi chr22.fa \
                   -bed <(awk '{if(NR > 1 && $2 != "WHG" && $2 != "LBK" && $2 != "Ust_ishim" && $3 == 22){print "chr"$3"\t"$4"\t"$5}}' Data1_archaicfragments.txt | sort -k1n -k2n -k3n | bedtools merge -i -) \
                   -fo chr22_masked.fa
```

Finally, the file [`chr22_masked.fa`](chr22_masked.fa) would have both repetitive regions and archaic-introgressed regions masked.  

[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="Der"></a>
## 4. Derived allele counts per mutation type

With the data produced, we can count the number of derived alleles accumulated after the out of Africa and outside repeatitive regions and outside archaic-introgression regions. These, will be distributed into the different mutation types taking into account the mutation context (5' and 3' base pairs surrounding the focal SNP).

[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="Ref"></a>
### 4.1. Load reference genomes

The human reference chromosome 22 with masked regions and the chimpanzee reference chromosome 22 (provided by SGDP) are load into memory.


```python
def ref_fasta(ref, chrom):
    '''
    This function will read a fasta file (ref), look for a chromosome (chrom) and output
    a string with the fasta sequence.
    '''
    fasta = ""
    with open("{}".format(ref), "r") as file:
        for line in file: 
            if line[0] == ">":
                c = line.strip().split(" ")[0][1:]
            else:
                if c == chrom or c == "chr"+chrom:
                    fasta = fasta + line.strip()
    return fasta

href  = ref_fasta("chr22_masked.fa", "22")
chimp = ref_fasta("Chimp.fa",        "22")
```


```python
href[16050000:16050010]
```




    'GATCTGATAA'



[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="Dat"></a>
### 4.2. ind, snp and geno files

From the `.ind` file, we extract three lists:

1. `ind` : with the individual names
2. `reg` : the regions the individuals with the same index in the `ind` list belong to
3. `afr` : index of the african individuals


```python
ind = []
sex = []
reg = []
afr = []

with open("SGDP_22.ind") as file:
    for i, line in enumerate(file):
        ind.append(line.strip().split()[0])
        sex.append(line.strip().split()[1])
        reg.append(line.strip().split()[2])
        if line.strip().split()[2] == "Africa":
            afr.append(i)
```


```python
print(ind[10:20]) #Subset of the ind list
print(reg[10:20]) #Corresponding subset of the reg list
print(afr[:10])   #First 10 african individual indexes
print(len(ind))   #Number of individuals
```

    ['S_Ami-2', 'S_Armenian-1', 'S_Armenian-2', 'S_Atayal-1', 'S_Balochi-1', 'S_Balochi-2', 'S_BantuHerero-1', 'S_BantuHerero-2', 'S_BantuKenya-1', 'S_BantuKenya-2']
    ['EastAsia', 'WestEurasia', 'WestEurasia', 'EastAsia', 'SouthAsia', 'SouthAsia', 'Africa', 'Africa', 'Africa', 'Africa']
    [16, 17, 18, 19, 20, 21, 30, 31, 52, 53]
    235


From the `.snp` file, we will construct a two-level dictionary with different information:

- The first level will be denoted by the index position of the SNP.
- The second level will have the following possible keys:
    1. "ref" : reference allele
    2. "alt" : alternative allele
    3. "pos" : genomic position of the SNP
    4. "anc" : ancestral allele (polarized with the chimp reference genome)
    5. "der" : derived allele (polarized with the chimp reference genome)
    6. "mut" : mutation type as fiv+anc+thr>der. E.g. ACG>G. Strand complentary mutation types are collapsed, being C and T the ancestral alleles. Eg. CGT>C is strand complentary to ACG>G, thus this mutation type will be denoted as ACG>G.

Another dictionary is constructed, `mut`, to keep which row indeces in the `.geno` file correspond to every mutation type.


```python
called_bases = ["A", "C", "T", "G"]
comp_base    = {"A" : "T",
                "T" : "A",
                "C" : "G",
                "G" : "C"}

def mut_type(fiv, anc, thr, der):
    '''
    Input:
        - fiv: 5' nucleotide
        - anc: ancestral nucleotide
        - thr: 3' nucleotide
        - der: derived nucleotide
    Output:
        - mutation type
    '''
    if anc in ["C", "T"]:
        return fiv+anc+thr+">"+der
    else:
        return comp_base[thr]+comp_base[anc]+comp_base[fiv]+">"+comp_base[der]
        

snp = defaultdict(lambda : defaultdict(lambda : None))
mut = defaultdict(lambda :  [])
        
with open("SGDP_22.snp") as file:
    for i, line in enumerate(file):
        _, _, _, s, r, a = line.strip().split()
        s = int(s)-1
        snp[i]["ref"] = r
        snp[i]["alt"] = a
        snp[i]["pos"] = s
        if href[s-1].upper() in called_bases and href[s] in called_bases and href[s+1].upper() in called_bases: #Checking if the position is not masked (Repetitive or archaic-introgressed region) and if the context (5' and 3' bases) are also called
            if chimp[s] == snp[i]["ref"]: #Polarizing the derived and ancestral allele in respect to the chimp allele
                snp[i]["anc"] = snp[i]["ref"]
                snp[i]["der"] = snp[i]["alt"]
            elif chimp[s] == snp[i]["alt"]: #Polarizing the derived and ancestral allele in respect to the chimp allele
                snp[i]["anc"] = snp[i]["alt"]
                snp[i]["der"] = snp[i]["ref"]
        if snp[i]["anc"]:
            snp[i]["mut"] = mut_type(href[s-1].upper(), snp[i]["anc"], href[s+1].upper(), snp[i]["der"]) #Define the mutation type
            mut[snp[i]["mut"]].append(i)
```


```python
print(mut["TCC>T"][:10]) # first 10 row-index corresponding to polymorphic sites which are TCC>T mutation types
print(snp[18])           #Dictionary entry example
print(len(snp))          #Number of SNPs
```

    [8791, 8793, 9010, 9021, 9077, 9091, 9129, 9157, 9172, 9216]
    defaultdict(<function <lambda>.<locals>.<lambda> at 0x2b24b284f820>, {'ref': 'A', 'alt': 'G', 'pos': 17053171, 'anc': 'A', 'der': 'G', 'mut': 'GTT>C'})
    345126


The geno file is imported as a numpy array


```python
geno = []

rev = {"0" : "2",
       "1" : "1",
       "2" : "0",
       "9" : "9"}

with open("SGDP_22.geno") as file:
    for i, line in enumerate(file):
        if snp[i]["anc"] == snp[i]["ref"]:
            geno.append([int(rev[g]) for g in line.strip()])
        else:
            geno.append([int(g) for g in line.strip()])

geno = np.array(geno, dtype = "int")
```


```python
print(geno)       #geno array
print(geno.shape) #geno array shape, if correct, equal to the number of SNPs and number of individuals
```

    [[2 2 2 ... 2 2 2]
     [2 1 1 ... 2 2 2]
     [2 2 2 ... 2 2 9]
     ...
     [2 2 9 ... 9 2 2]
     [2 2 9 ... 9 9 9]
     [2 9 9 ... 9 9 9]]
    (345126, 235)


[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="Fil"></a>
### 4.3. Filters

For every site, we applied certain filters:

1. Polymorphic in africa: in order to focus on mutations accumulated after the out of Africa, we filter out positions for which the derived allele is present in Africa. 
2. Callability per position: a polymorphic site for which less than 20% of the individuals are called, they are filtered out.
3. Masked, context noncallalbe or noncallable in chimp: if a position is in a masked region (repeat region or archaic-introgressed region), if its context (5' and 3' context) is noncallalbe or was not callalbe in chimps so it was not possible to polarize the alleles, then the position is going to be filtered out.

The genotypes of every individual for positions that don't pass all filters are going to be changed by the number 9 which encodes for noncallalble genotype.


```python
for i in range(geno.shape[0]):
    filter = False
    #1. CALLABILITY PER POSITION
    if np.sum(geno[i, 1:] == 9)/(geno.shape[1]-1) > 0.2: # more than 20% of the individuals called (disregarding the Href)
        filter = True
    
    #2. PRESENT IN AFRICA
    if np.sum(geno[i, afr][geno[i, afr] != 9]) > 0: # if there is presence of the derived allele in any of the african individuals
        filter = True

    #3. MASKED, CONTEXT NONCALLALBE, NONCALLABLE IN CHIMP
    if not snp[i]["anc"]: # if the snp does not have an associated ancestral allele, it means that one of the three conditions failed (check when the snp dictionary was constructed above)
        filter = True
        
    if filter:
        geno[i] = [9]*geno.shape[1] # if the site must be filtered out, the genotypes for every individual are changed to "9" (noncallable)

```


```python
geno
```




    array([[9, 9, 9, ..., 9, 9, 9],
           [9, 9, 9, ..., 9, 9, 9],
           [9, 9, 9, ..., 9, 9, 9],
           ...,
           [9, 9, 9, ..., 9, 9, 9],
           [9, 9, 9, ..., 9, 9, 9],
           [9, 9, 9, ..., 9, 9, 9]])



[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;

<a name="Mut"></a>
### 4.4. Mutation types per individual

The following bit of code outputs the counts per individual for every mutaition type.



```python
out = open("counts_22.txt", "w")

out.write("\t".join(["ind", "reg", "sex", "chrom", "fiv", "anc", "thr", "der", "count"])+"\n")

for fiv in called_bases:
    for anc in ["C", "T"]:
        for thr in called_bases:
            for der in called_bases:
                for i in range(len(ind)):
                    if reg[i] not in ["Africa", "hg19ref"]:
                        out.write("\t".join( [ind[i], reg[i], sex[i], "22", fiv, anc, thr, der, str(np.sum(geno[mut[fiv+anc+thr+">"+der], i][geno[mut[fiv+anc+thr+">"+der], i] != 9]))])+"\n")

out.close()
```


```bash
%%bash 

head counts_22.txt
```

    ind	reg	sex	chrom	fiv	anc	thr	der	count
    S_Abkhasian-1	WestEurasia	M	22	A	C	A	A	0
    S_Abkhasian-2	WestEurasia	M	22	A	C	A	A	1
    S_Adygei-1	WestEurasia	M	22	A	C	A	A	0
    S_Adygei-2	WestEurasia	F	22	A	C	A	A	0
    S_Albanian-1	WestEurasia	F	22	A	C	A	A	0
    S_Aleut-1	CentralAsiaSiberia	M	22	A	C	A	A	1
    S_Aleut-2	CentralAsiaSiberia	F	22	A	C	A	A	0
    S_Altaian-1	CentralAsiaSiberia	M	22	A	C	A	A	2
    S_Ami-1	EastAsia	M	22	A	C	A	A	0


By repeating this process for all chromosomes, we obtain `mutation_spectrum.txt`


```bash
%%bash 

head mutation_spectrum.txt
```

    ind	reg	chrom	fiv	anc	thr	der	count
    S_Abkhasian-1	WestEurasia	1	A	C	A	A	17
    S_Abkhasian-2	WestEurasia	1	A	C	A	A	13
    S_Adygei-1	WestEurasia	1	A	C	A	A	21
    S_Adygei-2	WestEurasia	1	A	C	A	A	24
    S_Albanian-1	WestEurasia	1	A	C	A	A	27
    S_Aleut-1	CentralAsiaSiberia	1	A	C	A	A	17
    S_Aleut-2	CentralAsiaSiberia	1	A	C	A	A	21
    S_Altaian-1	CentralAsiaSiberia	1	A	C	A	A	19
    S_Ami-1	EastAsia	1	A	C	A	A	16


This file was obtained by running [`workflow.py`](workflow.py) with the [gwf](https://docs.gwf.app/) (version 1.7.1) software, which is a workflow tool that takes care of job dependencies and keeping track of pipeline completion. 

<a name="Bib"></a>
## 5. References

1. Mallick, S., Li, H., Lipson, M. et al. The Simons Genome Diversity Project: 300 genomes from 142 diverse populations. Nature 538, 201–206 (2016). https://doi.org/10.1038/nature18964
2. https://github.com/DReichLab/cTools
3. Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, 15 March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033
4. https://docs.gwf.app/

[<img src="arrow.png" width="100" style="float: left;">](#H) &nbsp;

&nbsp;
