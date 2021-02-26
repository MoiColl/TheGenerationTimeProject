from gwf import AnonymousTarget, Workflow
gwf = Workflow()

#A. TEMPLATES

##1.
'''
Unzips the downloaded reference fasta sequence
'''
def unzip_ref(chrom):
	inputs  = ["chr{chrom}.fa.gz".format(chrom = chrom)]
	outputs = ["chr{chrom}.fa".format(chrom = chrom)]
	options = {
			'memory'  : '1g',
			'walltime': '00:10:00',
			'account' : 'simons'
	}
	spec = '''
	zcat chr{chrom}.fa.gz > chr{chrom}.fa.tmp

	mv chr{chrom}.fa.tmp chr{chrom}.fa

	'''.format(chrom = chrom)

	return inputs, outputs, options, spec

##2.
'''
Masks archaic-introgressed regions in the human reference chromosome fasta file
'''
def mask(chrom):
	inputs  = ["chr{chrom}.fa".format(chrom = chrom), "Data1_archaicfragments.txt"]
	outputs = ["chr{chrom}_masked.fa".format(chrom = chrom)]
	options = {
			'memory'  : '1g',
			'walltime': '00:10:00',
			'account' : 'simons'
	}
	spec = '''
	bedtools maskfasta \
	               -fi chr{chrom}.fa \
                   -bed <(awk '{{if(NR > 1 && $2 != "WHG" && $2 != "LBK" && $2 != "Ust_ishim" && $3 == {chrom}){{print "chr"$3"\t"$4"\t"$5}}}}' Data1_archaicfragments.txt \
                         | sort -k1n -k2n -k3n \
                         | bedtools merge -i -) \
                   -fo chr{chrom}_masked.tmp

    mv chr{chrom}_masked.tmp chr{chrom}_masked.fa

	'''.format(chrom = chrom)

	return inputs, outputs, options, spec


##3.
'''
Counts the number of derived alleles per individual found outside repetitive regions and outside 
archaic-introgressed regions
'''
def derived_alleles(chrom):
	inputs  = ["chr{chrom}_masked.fa".format(chrom = chrom), "Chimp.fa", 
	           "SGDP_{chrom}.ind".format(chrom = chrom), 
	           "SGDP_{chrom}.snp".format(chrom = chrom), 
	           "SGDP_{chrom}.geno".format(chrom = chrom)]
	outputs = ["counts_{chrom}.txt".format(chrom = chrom)]
	if chrom in [str(c) for c in range(1, 23)]:
		outputs.append("CGenrichment_{}.txt".format(chrom))
	options = {
			'memory'  : '20g',
			'walltime': '00:10:00',
			'account' : 'simons'
	}
	spec = '''
	python derived_alleles.py {chrom}

	mv counts_{chrom}.tmp counts_{chrom}.txt
	mv CGenrichment_{chrom}.tmp CGenrichment_{chrom}.txt
	'''.format(chrom = chrom)

	return inputs, outputs, options, spec


##4.
'''
Joins all mutation counts from derived_alleles() template
'''
def join_derived_alleles():
	inputs  = ["counts_{chrom}.txt".format(chrom = chrom) for chrom in [c for c in range(1, 23)]+["X", "Y"]]
	outputs = ["mutation_spectrum.txt"]
	options = {
			'memory'  : '1g',
			'walltime': '00:01:00',
			'account' : 'simons'
	}
	spec = '''
	echo "" | awk '{{print "ind\treg\tsex\tchrom\tfiv\tanc\tthr\tder\tcounts"}}' > mutation_spectrum.tmp
	for chrom in `seq 1 22` X Y; do awk '{{if(NR>1){{print}}}}' counts_$chrom.txt >> mutation_spectrum.tmp; done
	mv mutation_spectrum.tmp mutation_spectrum.txt

	'''

	return inputs, outputs, options, spec

##5.
'''
Joins all CGenrichment files from derived_alleles() template
'''
def join_CGenrichment():
	inputs  = ["CGenrichment_{chrom}.txt".format(chrom = chrom) for chrom in [c for c in range(1, 23)]]
	outputs = ["CGenrichment.txt"]
	options = {
			'memory'  : '1g',
			'walltime': '00:01:00',
			'account' : 'simons'
	}
	spec = '''
	echo "" | awk '{{print "ind\treg\tsex\tchrom\tnonEnrCG\tnonEnrnonCG\tEnrCG\tEnrnonCG"}}' > CGenrichment.tmp
	for chrom in `seq 1 22`; do awk '{{if(NR>1){{print}}}}' CGenrichment_$chrom.txt >> CGenrichment.tmp; done
	mv CGenrichment.tmp CGenrichment.txt
	'''

	return inputs, outputs, options, spec


