#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

"""
Input file from vcftools:

vcftools --vcf <vcf_file_here> --extract-FORMAT-info GT

The output file wiil be: out.GT.FORMAT

Then, submit the out.GT.FORMAT to Summary4VCF.py:

Summary4VCF.py out.GT.FORMAT

"""

import sys

try:
	outgt=sys.argv[1]
	f = open(outgt,'r')
except:
	print(__doc__)
	sys.exit(1)

#Take the info for the samples
samples={}
listSamples=[]
line = f.readline()
listSamples = line.split()[2:]
for y in listSamples:
	samples[y]=[]

line = f.readline()

#Take the info for the loci (sites)
sites={}
listSites=[]
while line:
	chr = line.split()[0]
	pos = line.split()[1]
	listSites.append(chr+'_'+pos)
	listGT = line.split()[2:]
	#Fill the sites dict
	sites[chr+'_'+pos] = listGT
	
	#Fill the samples dict
	a=0
	while a < len(listGT):
		samples[listSamples[a]].append(listGT[a])
		a+=1
	
	line = f.readline()

#Output the info in 2 files
#one for the sites
out1=open('Summary_By_Sites_python.txt','w')
out1.write('Chr_Pos\tTotal Number of Samples\tNumber of homozygotes\tNumber of heterozygotes\tNumber of other genotypes\tNumber of missing data\n')

for s in listSites:
	nb_homo1=sites[s].count('0/0')
	nb_homo2=sites[s].count('1/1')
	nb_het1=sites[s].count('1/0')
	nb_het2=sites[s].count('0/1')
	nb_miss=sites[s].count('./.')
	nb_others=len(sites[s])-nb_homo1-nb_homo2-nb_het1-nb_het2-nb_miss
#	out1.write(s+'\t'+str(len(sites[s]))+'\t'+str(nb_homo1+nb_homo2+nb_het1+nb_het2+nb_miss)+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(nb_miss)+'\n')
	out1.write(s+'\t'+str(len(sites[s]))+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(nb_others)+'\t'+str(nb_miss)+'\n')


#the other for the samples
out2=open('Summary_By_Samples_python.txt','w')
out2.write('Samples\tTotal Number of Sites\tNumber of homozygotes\tNumber of heterozygotes\tNumber of other genotypes\tNumber of missing data\n')
for s in listSamples:
	nb_homo1=samples[s].count('0/0')
	nb_homo2=samples[s].count('1/1')
	nb_het1=samples[s].count('1/0')
	nb_het2=samples[s].count('0/1')
	nb_miss=samples[s].count('./.')
	nb_others=len(samples[s])-nb_homo1-nb_homo2-nb_het1-nb_het2-nb_miss
#	out2.write(s+'\t'+str(len(samples[s]))+'\t'+str(nb_homo1+nb_homo2+nb_het1+nb_het2+nb_miss)+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(nb_miss)+'\n')
	out2.write(s+'\t'+str(len(samples[s]))+'\t'+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(nb_others)+'\t'+str(nb_miss)+'\n')

f.close()
out1.close()
out2.close()