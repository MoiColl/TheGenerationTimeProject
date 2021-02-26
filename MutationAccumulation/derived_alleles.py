#A. Importing
import sys
import numpy as np
from collections import defaultdict

#B. Global varaibles
called_bases = ["A", "C", "T", "G"]
comp_base    = {"A" : "T",
                "T" : "A",
                "C" : "G",
                "G" : "C",
                "X" : "X"}


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
            if chrom == "Y" and line.strip().split()[0] in ["S_Biaka-1", "S_Biaka-2", "S_Dinka-2", "S_Ju_hoan_North-1", "S_Ju_hoan_North-2", "S_Ju_hoan_North-3", "S_Mbuti-3"]:
                afr.append(i)
            if chrom != "Y" and line.strip().split()[2] == "Africa":
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
def inXDE(pos):
    '''
    Input:
        - pos   : locus genomic position
    Output:
        - boolean : If the locus is in X chromosome DEgenerate region (XDE), returns True.
    '''
    XDE = [[ 2649374,  2917723], [ 6616339,  7472224], 
           [14071703, 16095786], [16170060, 17986473], 
           [18016663, 18271273], [18537443, 19567356], 
           [21031901, 22216158], [22512750, 23497632]]

    for reg in XDE:
        start, end = reg
        if start <= pos and end > pos:
            return True
    return False

#C.5.
def get_snp(chrom, href, anc):
    '''
    By reading the .snp file and inputing the masked human reference fasta and the ancestral reference fasta
    get_snp() constructs two dictionaries:
    1. snp : this dictionary two-level dictionary stores for each snp the following information:
        1.1. i : SNP index. this is an integer that denotes the order of apearence of the SNPs in the .snp file
                 and consequently the index of the SNP in the .geno file
            1.1.1. ref : reference allele
            1.1.2. alt : alternative allele
            1.1.3. pos : genomic SNP position
            1.1.4. anc : ancestral allele (polarized according to chimp or 6 primates EPO ancestral call). If the 
                         SNP has no ancestral allele (snp[i]["anc"] == None) the SNP has faild for either:
                            - SNP in a repetitive region
                            - SNP in an archaic-introgressed region
                            - SNP with noncallalbe context sequence (5' and 3' contiguous nucleotides)
                            - SNP has noncallable homologous nucleotide in Chimp reference sequence
                            - SNP has homologous nucleotide in Chimp reference sequence different than ref and alt 
                              alleles
            1.1.5. der : derived allele (polarized according to chimp or 6 primates EPO ancestral call)
            1.1.6. mut : mutation type as as fiv+anc+thr>der (fiv and thr denote the 5' and 3' contiguous nucleotides 
                         of the SNP in question). E.g. ACG>G. Strand complentary mutation types are collapsed, 
                         being C and T the ancestral alleles. Eg. CGT>C is strand complentary to ACG>G, thus this 
                         mutation type will be denoted as ACG>G. For chromosome Y, mutation context is "X" since
                         we don't filter for context. Then, ACG>G mutation would be denoted as XCX>G.
    2. mut : dictionary of the 96 possible mutation types (denoted as fiv+anc+thr>der) each linked to a list with
             the indeces of the SNPs with the corresponding mutation type. For chromosome Y it's only 6 types.
    3. pos : numpy array with the genomic position of each snp

    4. CGp : dictionary which stores in a numpy array the genomic SNP positions of:
        4.1 C>G    : C>G mutations
        4.2 nonC>G : Mutations which are not C>G
    5. CGi : same as CGp but for indeces (i) instead of genomic positions
    '''
    pos = []
    snp = defaultdict(lambda : defaultdict(lambda : None))
    mut = defaultdict(lambda :  [])
    CGp = defaultdict(lambda :  [])
    CGi = defaultdict(lambda :  [])

    with open("SGDP_{}.snp".format(chrom)) as file:
        for i, line in enumerate(file):
            _, _, _, s, r, a = line.strip().split()
            s = int(s)-1
            snp[i]["ref"] = r
            snp[i]["alt"] = a
            snp[i]["pos"] = s
            pos.append(s)

            context = False
            if chrom != "Y":   #If it's not chromosome Y
                if href[s-1].upper() in called_bases and href[s] in called_bases and href[s+1].upper() in called_bases:#Checking if the 
                #position is not masked (Repetitive or archaic-introgressed region) and if the context (5' and 3' bases) are also called. 
                    context = True
            else:              #If it's chromosome Y
                if href[s] in called_bases and inXDE(s): #Check that it's a callable base pair and that it's in the XDE region.
                    context = True
            
            if context:
                if anc[s] == snp[i]["ref"]: #Polarizing the derived and ancestral allele in respect to the chimp allele or 6 primates EPO ancestral call
                    snp[i]["anc"] = snp[i]["ref"]
                    snp[i]["der"] = snp[i]["alt"]
                elif anc[s] == snp[i]["alt"]: #Polarizing the derived and ancestral allele in respect to the chimp allele or 6 primates EPO ancestral call
                    snp[i]["anc"] = snp[i]["alt"]
                    snp[i]["der"] = snp[i]["ref"]
                if snp[i]["anc"]:
                    if chrom != "Y":
                        snp[i]["mut"] = mut_type(href[s-1].upper(), snp[i]["anc"], href[s+1].upper(), snp[i]["der"]) #Define the mutation type
                    else:
                        snp[i]["mut"] = mut_type("X",               snp[i]["anc"], "X"              , snp[i]["der"])
                    
                    mut[snp[i]["mut"]].append(i)
                    if chrom not in ["Y", "X"]:
                        if snp[i]["mut"][1] == "C" and snp[i]["mut"][4] == "G":
                            CGp["C>G"].append(s)
                            CGi["C>G"].append(i)
                        else:
                            CGp["nonC>G"].append(s)
                            CGi["nonC>G"].append(i)


    pos           = np.array(pos)
    CGi["C>G"]    = np.array(CGi["C>G"])
    CGi["nonC>G"] = np.array(CGi["nonC>G"])
    CGp["C>G"]    = np.array(CGp["C>G"])
    CGp["nonC>G"] = np.array(CGp["nonC>G"])

    return pos, snp, mut, CGp, CGi


#C.6.
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

    chY = {"0" : "0",
           "1" : "9",
           "2" : "1",
           "9" : "9"}

    rcY = {"0" : "1",
           "1" : "9",
           "2" : "0",
           "9" : "9"}
    
    with open("SGDP_{}.geno".format(chrom)) as file:
        for i, line in enumerate(file):
            
            if chrom != "Y":
                if snp[i]["anc"] == snp[i]["ref"]:
                    geno.append([int(rev[g]) for g in line.strip()])
                else:
                    geno.append([int(g)      for g in line.strip()])

            else:
                if snp[i]["anc"] == snp[i]["ref"]:
                    geno.append([int(rcY[g]) for g in line.strip()])
                else:
                    geno.append([int(chY[g]) for g in line.strip()])
    
    geno = np.array(geno, dtype = "int")

    return geno

def get_ind_call(chrom, sex, ind):
    '''
    Returns the indeces of individuals which are going to be checked their callability for each loci
    to asess if a polymorphic site can be considered callable or not. This is beacuse, for instance, 
    Href is alway callable and has to be excluded; or for example, if we test a loci in chromosome Y
    we must also exclude all female indiviudals. 
    '''
    if chrom == "X":
        # For chromosome X, we exclude all males and the Href
        return np.where( (np.array(sex) == "F") * (np.array(ind) != "Href") )[0]
    elif chrom == "Y":
        # For chromosome Y, we exclude all females since they don't have Y chromosome and thus they
        # are going to be always uncallable. We also remove some male individuals which didn't yeald
        # any Y chromosome polymorphism despite them being males. We don't consider Href neither.
        return np.where( (np.array(sex) == "M") * np.array([True if x not in ["Href", "S_Finnish-2", "S_Finnish-3", "S_Palestinian-2", "S_Mansi-1", "S_Masai-2"] else False for x in ind]) )[0]
    else:
        # For autosomes, we exclude the Href sequence since it's always callable and doesn't give info
        # if a site has been correctly sequenced or not.
        return np.where(                           np.array(ind) != "Href"  )[0]


#C.7.
def filt_geno(geno, snp, afr, ind_call):
    '''
    This function changes all genotypes of all individuals to noncallalbe (9) of positions that:

        1. More than 20% of the individuals are noncallable 
        2. The derived allele is found in an african individual
        3. Filter conditions explained in the get_snp() function
    '''
    for i in range(geno.shape[0]):
        filter = False
        #1. CALLABILITY PER POSITION   

        if np.sum(geno[i, ind_call] == 9)/(ind_call.shape[0]) > 0.2: # more than 20% of the individuals called (disregarding the Href)
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


#C.8.
def CGenrichment(chrom, CGp, CGi, geno):
    '''
    Function which reads the Jonssons et al 2017 annotation of 1 Mb regions enriched for C>G mutations
    and counts the number of derived alleles that are C>G or nonC>G in C>G enriched in each window type.
    This is done for all individuals.
    '''
    CGcounts = defaultdict(lambda : defaultdict(lambda : np.zeros(len(ind))))
    with open("41586_2017_BFnature24018_MOESM3_ESM.txt") as file:
        for i, line in enumerate(file):
            if i:
                chrom_CG, window, CG_enriched = line.strip().split()
                window = int(window)
                if chrom == chrom_CG[3:]:
                    for muttype in ["C>G", "nonC>G"]:
                        pos_in_range = np.where((CGp[muttype] >= (window*1000000)) * (CGp[muttype] < ((window+1)*1000000)))[0]
                        CGcounts[CG_enriched][muttype] += np.sum(np.where(geno[CGi[muttype][pos_in_range], :] == 9, 0, geno[CGi[muttype][pos_in_range], :]), axis = 0) 
    return CGcounts

#D. Code

chrom                   = sys.argv[1]
href                    = ref_fasta("chr{}_masked.fa".format(chrom), chrom)
if chrom != "Y":
    anc                 = ref_fasta("Chimp.fa",                      chrom)
else:
    anc                 = ref_fasta("PanTro6inhg19coor_chrY.fa",    chrom)
ind, reg, sex, afr      = get_ind(chrom)
pos, snp, mut, CGp, CGi = get_snp(chrom, href, anc)
ind_call                = get_ind_call(chrom, sex, ind)
geno                    = filt_geno(get_geno(chrom, snp), snp, afr, ind_call)

#This section of the code outputs all the counts per individual, per chromosome and per mutation type

out = open("counts_{}.tmp".format(chrom), "w")
out.write("\t".join(["ind", "reg", "sex", "chrom", "fiv", "anc", "thr", "der", "counts"])+"\n")

if chrom == "Y":
    thr_fiv_called_bases = ["X"]
else:
    thr_fiv_called_bases = called_bases

for fiv in thr_fiv_called_bases:
    for anc in ["C", "T"]:
        for thr in thr_fiv_called_bases:
            for der in called_bases:
                if anc != der:
                    for i in range(len(ind)):
                        if reg[i] not in ["Africa", "hg19ref"]:
                            if (chrom not in ["X", "Y"]) or (chrom == "X" and sex[i] == "F") or (chrom == "Y" and sex[i] == "M" and ind[i] not in ["S_Finnish-2", "S_Finnish-3", "S_Palestinian-2", "S_Mansi-1", "S_Masai-2"]):
                                out.write("\t".join( [ind[i], reg[i], sex[i], chrom, fiv, anc, thr, der, str(np.sum(geno[mut[fiv+anc+thr+">"+der], i][geno[mut[fiv+anc+thr+">"+der], i] != 9]))])+"\n")

out.close()

#This section of the code outputs the counts of C>G and nonC>G mutations in C>G enriched and no C>G-enriched regions per individual. This is only done
#for autosomes.

if chrom not in ["X", "Y"]:

    CGcounts = CGenrichment(chrom, CGp, CGi, geno)

    out = open("CGenrichment_{}.txt".format(chrom), "w")
    out.write("\t".join(["ind", "reg", "sex", "chrom", "nonEnrC>G", "nonEnrnonC>G", "EnrnonC>G", "EnrnonC>G"])+"\n")

    for i in range(len(ind)):
        if reg[i] not in ["hg19ref", "Africa"]:
            out.write("\t".join([ind[i], reg[i], sex[i], chrom, str(int(CGcounts["0"]["C>G"][i])), 
                                                                str(int(CGcounts["0"]["nonC>G"][i])), 
                                                                str(int(CGcounts["1"]["C>G"][i])), 
                                                                str(int(CGcounts["1"]["nonC>G"][i]))])+"\n")
    out.close()



