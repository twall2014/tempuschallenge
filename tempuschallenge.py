import requests

def tempuschallenge():
	data_old = open('Challenge_data.vcf','r')
	data_new = open('Annotated_data.vcf','w')
	data = data_old.read().split('\n')[:-1]
	for line in data:
		if line[0] != '#':
			# split line into each possible element
			lineData = line.split('\t')

			# get type of variant
			typeOfVariant = lineData[-4].split('=')[-1]

			#parse flags of variant to get additional metadata
			flags = lineData[-2].split(':')
			depth = flags[2] # sequencing depth

			#determine most deleterious variant
			potentialVars = lineData[4].split(',')
			var = min(potentialVars,key=len)
			varIndex = potentialVars.index(var) #get index of deleterious variant for future 

			# get percentages of reference/alternate alleles
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
			response = requests.get('http://exac.hms.harvard.edu/rest/variant/'+'-'.join(varMetadata))
			response = response.text.encode('utf-8').split()
			try:
				queryIndex = response.index(key)
				allele_freq = response[queryIndex+1]
				allele_freq = float(allele_freq[:-1])
			except ValueError: #variant not found in ExAC, or no info found for allele frequency
				allele_freq = 0.0 #we assume allele hasn't been found yet

			newData = [typeOfVariant,depth,str(supportingVar),str(percentSupportingVar),str(percentSupportingRef),str(allele_freq)]

			line += '\t' + ':'.join(newData)
			data_new.write(line + '\n')
		else:
			if line[1] != '#':
				line += '\t' + 'alleleinfo'
			data_new.write(line + '\n')
			
if __name__ == '__main__':
	output_challenge()

