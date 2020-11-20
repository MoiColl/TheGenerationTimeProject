#A. Importing
import sys
import numpy as np
from collections import defaultdict

#B. Global varaibles
called_bases = ["A", "C", "T", "G"]
comp_base    = {"A" : "T",
                "T" : "A",
                "C" : "G",
                "G" : "C"}


#C. Functions

#C.1.
def ref_fasta(ref, chrom):
    '''
    This function will read a fasta file (ref), look for a chromosome (chrom) and output
    a string with the fasta sequence for the chromosome required
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

#C.2.
def get_ind(chrom):
    '''
    Reads the .ind file for the required chromosome and returns four lists:
    1. ind : list with the individual name/code
    2. sex : list with the individual sex (M = Male; F = Female)
    3. reg : list with the individual region (WestEurasia, SouthAsia, America, CentralAsiaSiberia, EastAsia, Africa)
    4. afr : list of the indeces in the previous lists corresponding to African individuals
    '''
    ind = []
    reg = []
    sex = []
    afr = []

    with open("SGDP_{}.ind".format(chrom)) as file:
        for i, line in enumerate(file):
            ind.append(line.strip().split()[0])
            sex.append(line.strip().split()[1])
            reg.append(line.strip().split()[2])
            if line.strip().split()[2] == "Africa":
                afr.append(i)

    return ind, reg, sex, afr

#C.3.
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

#C.4.
def get_snp(chrom, href, chimp):
    '''
    By reading the .snp file and inputing the masked human reference fasta and the chimpanzee reference fasta
    get_snp() constructs two dictionaries:
    1. snp : this dictionary two-level dictionary stores for each snp the following information:
        1.1. i : SNP index. this is an integer that denotes the order of apearence of the SNPs in the .snp file
                 and consequently the index of the SNP in the .geno file
            1.1.1. ref : reference allele
            1.1.2. alt : alternative allele
            1.1.3. pos : genomic SNP position
            1.1.4. anc : ancestral allele (polarized according to chimp). If the SNP has no ancestral allele
                         (snp[i]["anc"] == None) the SNP has faild for either:
                            - SNP in a repetitive region
                            - SNP in an archaic-introgressed region
                            - SNP with noncallalbe context sequence (5' and 3' contiguous nucleotides)
                            - SNP has noncallable homologous nucleotide in Chimp reference sequence
                            - SNP has homologous nucleotide in Chimp reference sequence different than ref and alt 
                              alleles
            1.1.5. der : derived allele (polarized according to chimp)
            1.1.6. mut : mutation type as as fiv+anc+thr>der (fiv and thr denote the 5' and 3' contiguous nucleotides 
                         of the SNP in question). E.g. ACG>G. Strand complentary mutation types are collapsed, 
                         being C and T the ancestral alleles. Eg. CGT>C is strand complentary to ACG>G, thus this 
                         mutation type will be denoted as ACG>G.
    2. mut : dictionary of the 96 possible mutation types (denoted as fiv+anc+thr>der) each linked to a list with
             the indeces of the SNPs with the corresponding mutation type
    '''
    snp = defaultdict(lambda : defaultdict(lambda : None))
    mut = defaultdict(lambda :  [])
            
    with open("SGDP_{}.snp".format(chrom)) as file:
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
    return snp, mut

#C.5.
def get_geno(chrom, snp):
    '''
    From the .geno file extracts the matrix of genotypes for all individuals for all SNPs. With the
    snp dictionary, it polarizes genotypes such that they are informative of the number of derived 
    allele copies for every SNP (reversed from the original notation of the .geno file).

    '''
    geno = []

    rev = {"0" : "2",
           "1" : "1",
           "2" : "0",
           "9" : "9"}
    
    with open("SGDP_{}.geno".format(chrom)) as file:
        for i, line in enumerate(file):
            if snp[i]["anc"] == snp[i]["ref"]:
                geno.append([int(rev[g]) for g in line.strip()])
            else:
                geno.append([int(g) for g in line.strip()])
    
    geno = np.array(geno, dtype = "int")

    return geno

def filt_geno(geno, snp, afr):
    '''
    This function changes all genotypes of all individuals to noncallalbe (9) of positions that:

        1. More than 20% of the individuals are noncallable 
        2. The derived allele is found in an african individual
        3. Filter conditions explained in the get_snp() function
    '''
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

    return geno


#D. Code
chrom              = sys.argv[1]
href               = ref_fasta("chr{}_masked.fa".format(chrom), chrom)
chimp              = ref_fasta("Chimp.fa",                      chrom)
ind, reg, sex, afr = get_ind(chrom)
snp, mut           = get_snp(chrom, href, chimp)
geno               = filt_geno(get_geno(chrom, snp), snp, afr)

#This section of the code outputs all the counts per individual, per chromosome and per mutation type
out = open("counts_{}.tmp".format(chrom), "w")
out.write("\t".join(["ind", "reg", "sex", "chrom", "fiv", "anc", "thr", "der", "counts"])+"\n")

for fiv in called_bases:
    for anc in ["C", "T"]:
        for thr in called_bases:
            for der in called_bases:
                for i in range(len(ind)):
                    if reg[i] not in ["Africa", "hg19ref"]:
                        out.write("\t".join( [ind[i], reg[i], sex[i], chrom, fiv, anc, thr, der, str(np.sum(geno[mut[fiv+anc+thr+">"+der], i][geno[mut[fiv+anc+thr+">"+der], i] != 9]))])+"\n")

out.close()





