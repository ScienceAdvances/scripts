import pandas as pd
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
	system("rm "+talign)
	return total,RiP,f'{score*100:.2f}%'

aligned_reads_bam = ["Result/BAM/KO_1.bam","Result/BAM/KO_2.bam","Result/BAM/WT_1.bam","Result/BAM/WT_2.bam"]
macs2_peak = ["Result/MACS/KO_1_peaks.narrowPeak","Result/MACS/KO_2_peaks.narrowPeak","Result/MACS/WT_1_peaks.narrowPeak","Result/MACS/WT_2_peaks.narrowPeak"]

total_reads_in_peaks=[]
total_mapped_reads=[]
frip=[]
for x,y in zip(aligned_reads_bam,macs2_peak):
    a1,a2,a3=FRiPscore(x,y)
    total_mapped_reads.append(a1)
    total_reads_in_peaks.append(a2)
    frip.append(a3)

df=pd.DataFrame(dict(total_reads_in_peaks=total_reads_in_peaks,total_mapped_reads=total_mapped_reads,FRiP=frip))
df.index = pd.Index(["KO","KO","WT","WT"],name='Sample')
df.to_csv('Report/table3.csv')
