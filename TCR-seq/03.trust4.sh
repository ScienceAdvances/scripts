perl BuildImgtAnnot.pl Homo_sapien > IMGT+C.fa
grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq > bcr_tcr_gene_name.txt  
perl BuildDatabaseFa.pl \
/home/alex/Alex/DataHub/GRCh38.primary_assembly.genome.fa \
/home/alex/Alex/DataHub/gencode.v48.primary_assembly.annotation.gtf \
bcr_tcr_gene_name.txt > hg38_bcrtcr.fa

sudo chmod u+x /home/alex/Alex/APP/TRUST4-master/run-trust4 

mamba env create -n fastp fastp -c bioconda -y
mamba activate fastp
ls Result/00.raw_data/*.gz | while read i
do
base=`basename $i`
name=${base%%.*}
echo $name
fastp \
    -i Result/00.raw_data/${name}.fastq.gz \
    -o Result/01.clean_data/${name}.fastq.gz \
    -w 16 -h /Result/fasatp/${name}.html \
    -j /Result/fasatp/${name}.json
done

ls Result/01.clean_data/*.gz | while read i
do
base=`basename $i`
name=${base%%.*}
echo $name
/home/alex/Alex/APP/TRUST4-master/run-trust4 \
-f src/hg38_bcrtcr.fa \
--ref src/human_IMGT+C.fa \
-u Result/01.clean_data/${name}.fastq.gz \
-o ${name} \
--od Result/03.TRUST4 \
-t 8
done
