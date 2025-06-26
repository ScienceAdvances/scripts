import argparse
import tempfile
import subprocess
from os import system
import uuid
import sys
def FRiPscore(bamFile, peaksFile, temporal=tempfile.gettempdir()):
	tempID = str(uuid.uuid4())
	score = None
	total = None
	#convert BAM to BED
	talign = temporal+"/"+tempID+".tagAlign"
	try:
		command = "bedtools bamtobed -i "+bamFile+" | awk 'BEGIN{OFS=\"\\t\"}{$4=\"N\";$5=\"1000\";print $0}' > "+talign
	except:
		print("Error: You need to install the bedtools suite. Exiting...")
		exit()
	system(command)
	#total reads in BAM
	command = "samtools view -c "+bamFile
	try:
		output = subprocess.check_output(command, shell=True)
		total = int(output)
	except:
		print("Error: You need to install the samtools suite. Exiting...")
		exit()
	#count reads in peak regions
	command = "bedtools sort -i "+peaksFile+" | bedtools merge -i stdin | bedtools intersect -u -a "+talign+" -b stdin | wc -l"
	RiP = int(subprocess.check_output(command, shell=True))
	if total > 0:
		score = RiP/total
	print(bamFile,"evaluated with",peaksFile,"has a FRiP score of: "+str(score))
	system("rm "+talign)

FRiPscore(args.bam,args.bed, args.temp)
