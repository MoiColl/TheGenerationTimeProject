#A. Importing
from templates import *
gwf = Workflow()


#B. Code
# for chrom in [str(c) for c in range(1:23)]+["X", "Y"]:
for chrom in [str(c) for c in range(1,23)]:
	gwf.target_from_template("unz_{}".format(chrom),
		unzip_ref(chrom = chrom))
	gwf.target_from_template("mas_{}".format(chrom),
		mask(chrom = chrom))
	gwf.target_from_template("der_{}".format(chrom),
		derived_alleles(chrom = chrom))

gwf.target_from_template("mut_spec",
	join_derived_alleles())
