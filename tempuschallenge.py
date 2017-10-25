import requests
from subprocess import call

def tempuschallenge():
	# This challenge takes in Challenge_data.vcf and returns a CSV file containing annotated information about the genome.
	# Output file is Annotated_data.csv.
	# Information in CSV in order:
	# - Chromosome
	# - Index
	# - Reference sequence
	# - Variant sequence
	# - Type of mutation
	# - Depth of sequence coverage
	# - Proportion of reads supporting variant allele
	# - Proportion of reads supporting reference allele
	# - ExAC allele frequency of variant

	# Run ANNOVAR on input data to get exact allele functions
	call(['perl', '../annovar/table_annovar.pl', 'Challenge_data.vcf', '../annovar/humandb/', '-buildver', 'hg19', '-out', 'annotated', '-remove', '-protocol', 'refGene,cytoBand,exac03', '-operation', 'g,r,f', '-nastring','.','-vcfinput'])
	
	# Open files for processing 
	data_old = open('annotated.hg19_multianno.vcf','r')
	data_new = open('Annotated_data.csv','w')
	data = data_old.read().split('\n')[:-1]

	# columns for new data table
	data_new.write("chrom,index,reference,variant,type,depth,percentSupportingVariant,percentSupportingReference,ExACAlleleFreq\n")
	for line in data:
		if line[0] != '#': # ignore input flags
			# split line into each possible element
			lineData = line.split('\t')

			#parse ending flags of variant to get additional metadata
			flags = lineData[-2].split(':')
			depth = flags[2] # sequencing depth

			#determine most deleterious variant
			potentialVars = lineData[4].split(',')
			var = min(potentialVars,key=len)
			varIndex = potentialVars.index(var) #get index of deleterious variant for future reference

			# get type of variant from ANNOVAR annotation
			variantFlags = lineData[-4].split(';')
			checkIndex = 0
			needExactFunction = False
			for f in variantFlags:
				if "Func.refGene" in f:
					if checkIndex == varIndex: # want to make sure correct variant info is returned when VCF indicates more than one possible variant
						if "ExonicFunc.refGene" in f and needExactFunction: # if mutation is exonic, need exact effect of mutation on coding region
							typeOfVariant = f.split('=')[1]
							if typeOfVariant == "stopgain": # clarifying ANNOVAR teminology
								typeOfVariant = "nonsense"
							elif typeOfVariant == "stoploss":
								typeOfVariant = "stop_codon_loss"
							elif typeOfVariant == "nonsynonymous_SNV":
								typeOfVariant = "missense"
							elif typeOfVariant == "synonymous_SNV":
								typeOfVariant = "silent"
							break
						else: 
							typeOfVariant = f.split('=')[1]
							if typeOfVariant == "exonic":
								needExactFunction = True
							else: 
								break
					else:
						checkIndex += 1
					
			# get percentages of reference/alternate alleles compared to read depth
			supportingVar = int(flags[4]);
			supportingRef = int(flags[6].split(',')[varIndex]);
			percentSupportingVar = float(supportingVar)/int(depth)
			percentSupportingRef = float(supportingRef)/int(depth)

			# metadata for ExAC query
			chrom = lineData[0]
			index = lineData[1]
			ref = lineData[3]
			varMetadata = [chrom,index,ref,var]

			# get relevant ExAC allele frequency info
			key = '"allele_freq":'
			response = requests.get('http://exac.hms.harvard.edu/rest/variant/'+'-'.join(varMetadata)) # query ExAC server for info
			response = response.text.encode('utf-8').split()
			try:
				queryIndex = response.index(key)
				allele_freq = response[queryIndex+1]
				allele_freq = float(allele_freq[:-1])
			except ValueError: #variant not found in ExAC, or no info found for allele frequency
				allele_freq = 0.0 #we assume allele hasn't been found yet

			# combine all data in line and write to file
			newData = varMetadata + [typeOfVariant,depth,str(supportingVar),str(percentSupportingVar),str(percentSupportingRef),str(allele_freq)]
			line = ','.join(newData)
			data_new.write(line + '\n')
	data_new.close()

if __name__ == '__main__':
	tempuschallenge()

