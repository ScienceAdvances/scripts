chmod u+x /home/alex/Alex/APP/vdjtools-1.2.1/vdjtools-1.2.1.jar
/home/alex/Alex/APP/vdjtools-1.2.1/vdjtools --help
install.packages("https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz")

brew install openjdk@8
export PATH="/home/linuxbrew/.linuxbrew/opt/openjdk@8/bin:$PATH"
export CPPFLAGS="-I/home/linuxbrew/.linuxbrew/opt/openjdk@8/include"


ls Result/TRUST4/*_report.tsv | while read i
do
base=`basename $i`
grep -v out_of_frame $i > Result/VDJtools/TRUST4/$base
done

vdjtools() {
    java -jar -Xmx15g /home/alex/Alex/APP/vdjtools-1.2.1/vdjtools-1.2.1.jar "$@"
}

vdjtools CalcBasicStats --metadata  Result/VDJtools/TRUST4_Metadata.tsv Result/VDJtools/CalcBasicStats/CalcBasicStats
vdjtools RarefactionPlot --plot --metadata Result/VDJtools/TRUST4_Metadata.tsv Result/VDJtools/RarefactionPlot/RarefactionPlot
awk 'NR>1 {print $3}' Result/VDJtools/CalcBasicStats/CalcBasicStats.basicstats.txt > Result/VDJtools/reads_list.txt
min_reads=$(sort -n Result/VDJtools/reads_list.txt | head -1)

vdjtools Correct --metadata Result/VDJtools/TRUST4_Metadata.tsv Result/VDJtools/Correct/Correct
vdjtools Decontaminate --metadata Result/VDJtools/Correct/TRUST4_Metadata.tsv Result/VDJtools/Decontaminate/Decontaminate
vdjtools DownSample -size $min_reads --metadata Result/VDJtools/Decontaminate/TRUST4_Metadata.tsv Result/VDJtools/DownSample/DownSample

metadata=Result/VDJtools/FilterNonFunctional/Metadata.tsv
vdjtools FilterNonFunctional -m Result/VDJtools/DownSample/TRUST4_Metadata.tsv -c Result/VDJtools/FilterNonFunctional/FilterNonFunctional
# Join samples into a single clonotype abundance matrix
vdjtools JoinSamples --plot --metadata $metadata Result/VDJtools/JoinSamples/JoinSamples
# Pool samples together
vdjtools PoolSamples --metadata $metadata Result/VDJtools/PoolSamples/PoolSamples


vdjtools CalcSegmentUsage --plot --metadata $metadata Result/VDJtools/CalcSegmentUsage/CalcSegmentUsage
vdjtools CalcDiversityStats --metadata $metadata Result/VDJtools/CalcDiversityStats/CalcDiversityStats
vdjtools CalcPairwiseDistances --plot --metadata $metadata Result/VDJtools/CalcPairwiseDistances/CalcPairwiseDistances
vdjtools CalcSpectratype --metadata $metadata Result/VDJtools/CalcSpectratype/CalcSpectratype
vdjtools CalcSpectratype --amino-acid --metadata $metadata Result/VDJtools/CalcSpectratype/CalcSpectratype

 
  