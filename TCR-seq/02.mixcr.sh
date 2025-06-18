mixcr=src/mixcr
sudo chmod -R u+x src/mixcr

ls Result/00.raw_data/*.gz | while read i
do
base=`basename $i`
name=${base%%.*}
echo $name
$mixcr analyze generic-amplicon \
    --species hsa \
    --rna \
    Result/00.raw_data/${name}.fastq.gz \
    Result/02.mixcr/${name} \
    --threads 16 \
    -Xmx15g \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary J \
    --assemble-clonotypes-by CDR3
done


mkdir -p Result/02.mixcr/clnx/ Result/02.mixcr/postanalysis
mkdir 
cp Result/02.mixcr/*clns Result/02.mixcr/clnx/

$mixcr postanalysis individual \
-Xmx15g \
--default-downsampling 'count-reads-auto' \
--default-weight-function read \
--chains TRB \
--tables Result/02.mixcr/postanalysis/tables.tsv \
--metadata meta.tsv \
--verbose \

$mixcr -Xmx15g exportPlots diversity \
    --plot-type boxplot \
    --primary-group Group \
    --primary-group-values "HC,Pre,Post" \
    Result/02.mixcr/postanalysis/result.json \
    diversity_facets.pdf
