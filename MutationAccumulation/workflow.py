#A. IMPORTING
from gwf import Workflow
from collections import defaultdict
import subprocess
import sys
import os
import yaml
import time
gwf = Workflow()

#B. TEMPLATES AND FUNCTIONS

##B.1. TEMPLATES

##B.1.1
'''
For the chromosome indicated, transforms the PhyloP files in wig format to bed format for latter use of them to mask different regions
'''
def phylop_wig2bed(chrom, completed):
	inputs  = []
	outputs = [completed]
	options = {
			'memory': '1g',
			'walltime': '00:40:00',
			'account': 'simons'
	}
	spec = '''

	echo "JOBID:" $PBS_JOBID

	python phylop_wig2bed.py {chrom}

	echo "Completed at "$(date) > {completed}'''.format(chrom = chrom, completed = completed)

	return inputs, outputs, options, spec


##B.1.2
'''
For the chromosome indicated, transforms the lower case leters (repeats) in fasta files to bed format for latter use of them to mask different regions
'''
def repeatmasker_fasta2bed(chrom, completed):
	inputs  = []
	outputs = [completed]
	options = {
			'memory': '1g',
			'walltime': '00:40:00',
			'account': 'simons'
	}
	spec = '''

	echo "JOBID:" $PBS_JOBID

	python repeatmasker_fasta2bed.py {chrom}

	echo "Completed at "$(date) > {completed}'''.format(chrom = chrom, completed = completed)

	return inputs, outputs, options, spec


##B.1.3
'''
For the chromosome indicated, get the joined neanderthal regions that all the individuals used in the pipeline have
'''
def arch_regions(prefix, chrom, completed):
	inputs  = []
	outputs = [completed]
	options = {
			'memory': '1g',
			'walltime': '00:10:00',
			'account': 'simons'
	}
	spec = """

	echo "JOBID:" $PBS_JOBID

	awk '{{if($2 == {chrom} && $16 == "Simons"){{print}}}}' /home/moicoll/GenerationTime/00Rawdata/ArchaicSegments.txt \
	| egrep -f <(awk '{{print $1}}' /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}.ind) \
	| awk '{{if($4>$3){{print "chr"$2"\t"$3"\t"$4}}}}' \
	| sort -k2,2n -k3,3n \
	| /com/extra/bedtools/2.25.0/bin/bedtools merge -i - > /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_archaic_chr{chrom}.bed

	echo "Completed at "$(date) > {completed}""".format(prefix = prefix, chrom = chrom, completed = completed)

	return inputs, outputs, options, spec


##B.1.4
'''
Sort all the regions to mask in the multiple files and create a single bedfile
'''
def masked_regions(phylop, repeatmask, arch_filt, prefix, chrom, prev_completed, completed):
	inputs  = prev_completed
	outputs = [completed]
	options = {
			'memory': '20g',
			'walltime': '01:00:00',
			'account': 'simons'
	}
	spec = """

	echo "JOBID:" $PBS_JOBID

	rm -f /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mask_chr{chrom}.bed

	if [ True == {repeatmask} ];
	then
		cat /home/moicoll/GenerationTime/faststorage/02MutationProfile/files/RepeatMasker/chr{chrom}.bed >> /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mask_chr{chrom}.bed
	fi;
	if [ True == {phylop} ];
	then
		cat /home/moicoll/GenerationTime/faststorage/02MutationProfile/files/phyloP100way/chr{chrom}.phyloP100way.bed >> /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mask_chr{chrom}.bed
	fi;
	if [ True == {arch_filt} ];
	then
		cat /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_archaic_chr{chrom}.bed >> /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mask_chr{chrom}.bed
	fi;

	sort -k2,2n -k3,3n /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mask_chr{chrom}.bed \
	| /com/extra/bedtools/2.25.0/bin/bedtools merge -i - > /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_repeatmasker_phylop_archaic_chr{chrom}.bed

	rm -f /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mask_chr{chrom}.bed

	echo "Completed at "$(date) > {completed}""".format(phylop = phylop, repeatmask = repeatmask, arch_filt = arch_filt, prefix = prefix, chrom = chrom, prev_completed = prev_completed, completed = completed)

	return inputs, outputs, options, spec


##B.1.5
'''
When a process is paralelized with multiple jobs, I use this template to flag that all the jobs for that specific task had been ran
'''
def all_completed(prev_completed, completed):
	inputs  = prev_completed
	outputs = [completed]
	options = {
			'memory': '1g',
			'walltime': '00:01:00',
			'account': 'simons'
	}
	spec = """
	
	echo "JOBID:" $PBS_JOBID


	echo 'Completed at '$(date) > {completed}""".format(prev_completed = prev_completed, completed = completed)

	return inputs, outputs, options, spec

##B.1.6
'''
Run the cpoly with the parameter and ind files created for each chromosome
'''
def cpoly(prefix, chrom, completed):
		inputs = ["/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}.par".format(prefix = prefix, chrom = chrom),
				  "/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}.ind".format(prefix = prefix)]
		outputs = [completed]
		options = {
				'memory': '2g',
				'walltime': '08:00:00',
				'account': 'simons'
		}
		spec = '''
		echo "JOBID:" $PBS_JOBID

		cpoly -p /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}.par

		sleep 2
		
		echo "Completed at "$(date) > {completed}'''.format(prefix = prefix, chrom = chrom, completed = completed)

		return inputs, outputs, options, spec

##B.1.6
'''
Run the cpoly with the parameter and ind files created for each chromosome
'''
# def cpoly_moi(prefix, chrom, start, end, qual, completed):
# 		inputs = ["/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}.par".format(prefix = prefix, chrom = chrom),
# 				  "/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}.ind".format(prefix = prefix)]
# 		outputs = [completed]
# 		options = {
# 				'memory': '8g',
# 				'walltime': '03:00:00',
# 				'account': 'simons'
# 		}
# 		spec = '''
# 		echo "JOBID:" $PBS_JOBID

# 		python extract_derived_alleles.py {prefix} {chrom} {start} {end} {qual}

# 		sleep 2
		
# 		echo "Completed at "$(date) > {completed}'''.format(prefix = prefix, chrom = chrom, start = start, end = end, qual = qual, completed = completed)

# 		return inputs, outputs, options, spec

##B.1.7
'''
Filter the raw data and obtain the 96 mutation spectrum
'''
def mutation_spectrum(chrom, prefix, mask, call, outgroup, not_pop_singleton, compare_freqs, comparison, freq, not_poly_africa, cluster, prev_completed, completed):
		inputs  = prev_completed
		outputs = [completed]
		options = {
				'memory': '36g',
				'walltime': '03:00:00',
				'account': 'simons'
		}
		spec = '''
		echo "JOBID:" $PBS_JOBID

		python mutation_spectrum.py {chrom} {prefix} {mask} {call} {outgroup} {not_pop_singleton} {compare_freqs} "{comparison}" {freq} {not_poly_africa} {cluster}

		sleep 2
		
		echo "Completed at "$(date) > {completed}'''.format(chrom = chrom, prefix = prefix, mask = mask, 
			                                                call = call, outgroup = outgroup, not_pop_singleton = not_pop_singleton, 
			                                                compare_freqs = compare_freqs, comparison = comparison, freq = freq,
			                                                not_poly_africa = not_poly_africa, cluster = cluster, completed = completed)

		return inputs, outputs, options, spec

##B.1.8
'''
Join all the mutation spectrum files
'''
def all_mutation_spectrum(prefix, prev_completed, completed):
	inputs  = prev_completed
	outputs = [completed]
	options = {
			'memory': '1g',
			'walltime': '00:10:00',
			'account': 'simons'
	}
	spec = """
	
	echo "JOBID:" $PBS_JOBID
 
	echo "" | awk '{{print "ind\tpop\tchrom\tfiv\tanc\tthr\tder\tcounts"}}' > /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mut_spectrum.txt

	cat /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_*_mut_spectrum.txt \
	| sort -k 1,1n -k 2,2n -k 3,3V >> /home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_mut_spectrum.txt

	echo 'Completed at '$(date) > {completed}""".format(prefix = prefix, completed = completed)

	return inputs, outputs, options, spec	


##B.1.9
'''
Mask fasta files for archaic samples to use as a single fasta file instead of both the hetfa + mask files since they give errors with Cpoly
'''
def ancient_archaic_samples_mask_fasta(sample, completed):
	inputs  = []
	outputs = [completed]
	options = {
			'memory': '1g',
			'walltime': '03:00:00',
			'account': 'simons'
	}
	spec = """
	
	echo "JOBID:" $PBS_JOBID

	python ancient_archaic_samples_mask_fasta.py {sample}

	echo 'Completed at '$(date) > {completed}""".format(sample = sample, completed = completed)

	return inputs, outputs, options, spec	



##B.2. FUNCTIONS

def chrom_list(prefix):
	chroms = []
	for chrom in range(1, 23):
		chroms.append(str(chrom))
	if prefix in ["SGDP_2", "X"]:
		chroms.append(str("X"))
	# chroms.append(str("Y"))
	return chroms

def parsing_variables_from_parameter_file(line, line_number):

	population_datasets = {"SGDP"  : ["WestEurasia", "EastAsia", "CentralAsiaSiberia", "SouthAsia", "America", "Oceania", "Africa", "Altai", "Denisova", "WHG", "LBK", "Ust_Ishim"],
						   "1000G" : ["WestEurasia", "EastAsia", "CentralAsiaSiberia", "SouthAsia", "America", "Oceania", "Africa"]} ##########################################

	prefix, dataset, phylop, repeatmask, arch_filt, arch_filt_pops, qual, call, not_pop_singleton, pops, pops_file, outgroup, compare_freqs, comparison, freq, not_poly_africa, cluster = line.strip().split(";")
	if dataset not in ["SGDP", "1000G"]:
		sys.exit("The 'dataset' is not any of the following options: 'SGDP', '1000G' in line {}".format(line_number))
	if phylop == "TRUE":
		phylop = True
	elif phylop == "FALSE":
		phylop = False
	else:
		sys.exit("I can't convert the 'phylop' filter parameter in line {}".format(line_number))
	if repeatmask == "TRUE":
		repeatmask = True
	elif repeatmask == "FALSE":
		repeatmask = False
	else:
		sys.exit("I can't convert the 'repeatmask' filter parameter in line {}".format(line_number))
	if arch_filt == "TRUE":
		arch_filt = True
		for pop in arch_filt_pops.split(","):
			if pop not in population_datasets[dataset]:
				sys.exit("In the 'arch_filt_pops' there is a population not registered in this dataset in line {}".format(line_number))
		arch_filt_pops = arch_filt_pops.split(",")
	elif arch_filt == "FALSE":
		arch_filt = False
	else:
		sys.exit("I can't convert the 'arch_filt' filter parameter in line {}".format(line_number))
	if dataset == "SGDP":
		try:
			qual = int(qual)
		except:
			sys.exit("I can't convert the 'qual' filter into a integer since {} dataset is used parameter in line {}".format(dataset, line_number))
	else:
		qual = "PASS"
	try:
		call = float(call)
		if call < 0 and call > 1:
			sys.exit("The callability percentage is not between 0 and 1 in line {}".format(line_number))
	except:
		sys.exit("The callability percentage is not a number between 0 and 1 in line {}".format(line_number))
	if not_pop_singleton == "TRUE":
		not_pop_singleton = True
	elif not_pop_singleton == "FALSE":
		not_pop_singleton = False
	else:
		sys.exit("I can't convert the 'not_pop_singleton' filter parameter in line {}".format(line_number))
	for pop in pops.split(","):
		if pop not in population_datasets[dataset]:
			sys.exit("In the 'arch_filt_pops' there is a population not registered in this dataset in line {}".format(line_number))
	pops = pops.split(",")
	if outgroup not in ["CHIMP", "ANC"]:
		sys.exit("In the 'outgroup' there is an outgroup that is not either CHIMP or ANC in line {}".format(line_number))
	if compare_freqs == "TRUE":
		compare_freqs = True
		if comparison not in [">", ">=", "<", "<="]:
			sys.exit("The operator in the field corresponding to the 'comparison' is not any of the > >= < <= operators in line {}".format(line_number))
		try:
			freq = float(freq)
			if freq < 0 or freq > 1:
				sys.exit("The variable 'freq' is not a numeric variable between 0 and 1 {}".format(line_number))
		except:
			sys.exit("The variable 'freq' is not a numeric variable between 0 and 1 {}".format(line_number))
	elif compare_freqs == "FALSE":
		compare_freqs = False
	else:
		sys.exit("I can't convert the 'compare_freqs' filter parameter in line {}".format(line_number))
	if not_poly_africa == "TRUE":
		not_poly_africa = True
	elif not_poly_africa == "FALSE":
		not_poly_africa = False
	else:
		sys.exit("I can't convert the 'not_poly_africa' filter parameter in line {}".format(line_number))
	if cluster == "FALSE":
		cluster = False
	else:
		try:
			cluster = int(cluster)
			if cluster <= 0:
				sys.exit("I can't convert the 'cluster' to a positive integer {}".format(line_number))
		except:
			sys.exit("I can't convert the 'cluster' to a positive integer or to the False boolean {}".format(line_number))

	return prefix, dataset, phylop, repeatmask, arch_filt, arch_filt_pops, qual, call, not_pop_singleton, pops, pops_file, outgroup, compare_freqs, comparison, freq, not_poly_africa, cluster

def create_par_file(prefix, chrom, qual):
	out = open("/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}.par".format(prefix = prefix, chrom = chrom), "w")
	out.write("indivname:\t/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}.ind\n".format(prefix = prefix))
	out.write("indivoutname:\t/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}_out.ind\n".format(prefix = prefix, chrom = chrom))
	out.write("snpoutname:\t/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}.snp\n".format(prefix = prefix,  chrom = chrom))
	out.write("genooutname:\t/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/{prefix}/{prefix}_{chrom}.geno\n".format(prefix = prefix, chrom = chrom))
	out.write("dbhetfa:\t/home/moicoll/GenerationTime/faststorage/02MutationProfile/files/hetfa.dblist\n")
	out.write("dbmask:\t/home/moicoll/GenerationTime/faststorage/02MutationProfile/files/mask.dblist\n")
	out.write("chrom:\t{chrom}\n".format(chrom = chrom))
	out.write("minfilterval:\t{qual}\n".format(qual = qual))
	out.write("polarize:\tHref\n")
	if chrom == "18":
		out.write("hipos:\t78000000\n")
	out.write("maxchrom:\t23\n")
	out.close()            
	

def create_ind_file(prefix):
	pop_dict = yaml.load(open("/home/moicoll/GenerationTime/faststorage/02MutationProfile/files/pop_yaml/"+prefix+".yaml"))
	ind_used = defaultdict(lambda : 0)
	for pop in pop_dict:
		if pop != "not_included":
			for ind in pop_dict[pop]:
				ind_used[ind] = 1
	ind_file = open("/home/moicoll/GenerationTime/faststorage/02MutationProfile/out/"+prefix+"/"+prefix+".ind", "w")
	with open("/home/moicoll/GenerationTime/faststorage/02MutationProfile/files/ind_pop_region_sex_correctsex.txt", "r") as file:
		for line in file:
			ind, subpop, region, sex = line.strip().split()
			if ind_used[ind]:
				ind_file.write("{}\t{}\t{}\n".format(ind, sex, region))
	ind_file.close()


#C. CODE
##C.1 VARIABLES
pwd                   = "/home/moicoll/GenerationTime/faststorage/02MutationProfile/"
files_dir             = pwd+"files/"
out_dir               = pwd+"out/" 
scripts_dir			  = pwd+"scripts/"
completed_dir         = out_dir+"COMPLETED/"
parameter_file        = pwd+"scripts/parameters.csv"
# callable_fraction_dir = pwd+"files/callable_fraction_per_ind/" #################################################################################################
# fasta_dir             = "/home/moicoll/simons/faststorage/data/cteam_lite_public3/FullyPublic/" ################################################################
os.system("mkdir -p {}".format(completed_dir))
# os.system("mkdir -p {}".format(callable_fraction_dir)) #########################################################################################################

# chr_len = {
# "1"  : 249250621, "2"  : 243199373, "3"  : 198022430, "4"  : 191154276, "5"  : 180915260,
# "6"  : 171115067, "7"  : 159138663, "8"  : 146364022, "9"  : 141213431, "10" : 135534747,
# "11" : 135006516, "12" : 133851895, "13" : 115169878, "14" : 107349540, "15" : 102531392,
# "16" : 90354753,  "17" : 81195210,  "18" : 78077248,  "19" : 59128983,  "20" : 63025520,
# "21" : 48129895,  "22" : 51304566,  "X"  : 155270560, "Y"  : 59373566}


##C.2 PREPARING ALL THE FILES NECESSARI
###C.2.1 Phylop files changing format to bed
phylop_wig2bed_completed = []
for chrom in chrom_list("X"):
	completed = completed_dir+"PHYLOP_WIG2BED_"+chrom+".COMPLETED"
	gwf.target_from_template("phylop_wig2bed_"+chrom, 
		phylop_wig2bed(chrom = chrom, completed = completed))
	phylop_wig2bed_completed.append(completed)
completed = completed_dir+"PHYLOP_WIG2BED.COMPLETED"
gwf.target_from_template("phylop_wig2bed", 
	all_completed(prev_completed = phylop_wig2bed_completed, completed = completed))

###C.2.2 Repeatmasker Fasta files changing format to bed
repeatmasker_fasta2bed_completed = []
for chrom in chrom_list("X"):
	completed = completed_dir+"REPEATMASKER_FASTA2BED_"+chrom+".COMPLETED"
	gwf.target_from_template("repeatmasker_fasta2bed_"+chrom, 
		repeatmasker_fasta2bed(chrom = chrom, completed = completed))
	repeatmasker_fasta2bed_completed.append(completed)
completed = completed_dir+"REPEATMASKER_FASTA2BED.COMPLETED"
gwf.target_from_template("repeatmasker_fasta2bed", 
	all_completed(prev_completed = repeatmasker_fasta2bed_completed, completed = completed))


###C.2.3 Masking fasta files for ancient and archaic samples
for sample in ["Altai", "Denisova", "Ust_Ishim", "Stuttgart", "Loschbour"]:
	completed = completed_dir+"ANCIENT_ARCHAIC_"+sample+"_MASK.COMPLETED"
	gwf.target_from_template("anci_arch_"+sample+"_mask", 
		ancient_archaic_samples_mask_fasta(sample = sample, completed = completed))

##C.3 RUNNING THE PIPELINE WITH THE MULTIPLE PARAMETERS
with open(parameter_file, "r") as file:
	for i, line in enumerate(file):
		if line[0] != "#":
			###C.3.1 GETTING VARIABLES FROM PARAMETER FILE
			prefix, dataset, phylop, repeatmask, arch_filt, arch_filt_pops, qual, call, not_pop_singleton, pops, pops_file, outgroup, compare_freqs, comparison, freq, not_poly_africa, cluster = parsing_variables_from_parameter_file(line, i)
			out_prefix_dir = out_dir+prefix+"/"
			os.system("mkdir -p {}".format(out_prefix_dir))
			completed_to_filter_step = []
			completed_cpoly = []

			###C.3.2 DECIDING WHICH DATASET ARE WE USING
			if dataset == "SGDP":
				####C.3.2.1 Simons Diversity Genome Project (SGDP)
				#####C.3.2.1.1 Preparing ind file to run cpoly 
				all_mutation_spectrum_completed = []
				create_ind_file(prefix)
				for chrom in chrom_list(prefix):
					#####C.3.2.1.2 Preparing par file to run cpoly 
					create_par_file(prefix, chrom, qual)
					#####C.3.2.1.3 Run cpoly
					completed = completed_dir+"CPOLY_"+chrom+"_"+prefix+".COMPLETED"
					if not os.path.isfile(completed):
						gwf.target_from_template("cpoly_"+chrom+"_"+prefix, 
							cpoly(prefix = prefix, chrom = chrom, completed = completed))
					completed_to_filter_step.append(completed)
					completed_cpoly.append(completed)
					#####C.3.2.1.4 Preparing the masking files
					mask = False
					completed_masking = []
					if phylop:
						completed_masking.append(completed_dir+"PHYLOP_WIG2BED_"+chrom+".COMPLETED")
						mask = True
					if repeatmask:
						completed_masking.append(completed_dir+"REPEATMASKER_FASTA2BED_"+chrom+".COMPLETED")
						mask = True
					if arch_filt:
						completed = completed_dir+"ARCH_REGIONS"+chrom+"_"+prefix+".COMPLETED"
						gwf.target_from_template("arch_regions_"+chrom+"_"+prefix, 
							arch_regions(prefix = prefix, chrom = chrom, completed = completed))
						completed_masking.append(completed)
						mask = True
					if mask:
						completed = completed_dir+"MASKED_REGIONS"+chrom+"_"+prefix+".COMPLETED"
						gwf.target_from_template("masked_regions_"+chrom+"_"+prefix, 
							masked_regions(phylop = str(phylop), repeatmask = str(repeatmask), arch_filt = str(arch_filt), prefix = prefix, chrom = chrom, prev_completed = completed_masking, completed = completed))
						completed_to_filter_step.append(completed)


				#####C.3.2.1.5 Mutation Profile by chromosome
				cpoly_completed_boolean = True
				for completed in completed_cpoly:
					if not os.path.isfile(completed):
						cpoly_completed_boolean = False
				if cpoly_completed_boolean:
					for chrom in chrom_list(prefix):
						completed = completed_dir+"MUT_SPECTRUM_"+chrom+"_"+prefix+".COMPLETED"
						gwf.target_from_template("mut_prof_"+chrom+"_"+prefix, 
							mutation_spectrum(chrom = chrom, prefix = prefix, mask = mask, call = call, outgroup = outgroup, not_pop_singleton = not_pop_singleton,
								compare_freqs = compare_freqs, comparison = comparison, freq = freq, not_poly_africa = not_poly_africa, cluster = cluster,
								prev_completed = completed_to_filter_step, completed = completed))
						all_mutation_spectrum_completed.append(completed)

					completed = completed_dir+"ALL_MUT_SPECTRUM_"+prefix+".COMPLETED"
					gwf.target_from_template("mut_prof_"+prefix, 
						all_mutation_spectrum(prefix = prefix, prev_completed = all_mutation_spectrum_completed, completed = completed))




			


