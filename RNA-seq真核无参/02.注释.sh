# cd Code/src
# cp -r resources $HOME/miniforge3/envs
# mkdir -p $HOME/miniforge3/envs/trinity/bin/util
# chmod -R u+x *.pl admin trinotateSeqLoader
# cp -r *.pl admin/ admin/util/*pl trinotateSeqLoader/ $HOME/miniforge3/envs/trinity/bin/util/
# cd ../../

# ======== Step0 构建数据库  ========
mkdir TrinotateDB; cd TrinotateDB
Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
cp ../uniprot_sprot.pep .
makeblastdb -in uniprot_sprot.pep -dbtype prot
diamond makedb --db uniprot_sprot.pep --in uniprot_sprot.pep
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
Trinotate --db Trinotate.sqlite --create --trinotate_data_dir .
cd ..

transcripts=transcripts.fa
proteins=transcripts.fa.transdecoder_dir/transcripts.fa.transdecoder.pep
uniprot_sprot_db=TrinotateDB/uniprot_sprot.dmnd
Pfam=TrinotateDB/Pfam-A.hmm
sqlitedb=TrinotateDB/Trinotate.sqlite

# ======== Step1 Run Transdecoder  ========
# （1）初步筛选ORF
TransDecoder.LongOrfs -S -t $transcripts
# （2）与已知蛋白数据库比对（可选）
diamond blastp \
    --query transcripts.fa.transdecoder_dir/longest_orfs.pep \
    --db $uniprot_sprot_db  \
    --max-target-seqs 1 \
    --outfmt 6 \
    --evalue 1e-5 \
    --threads 20  \
    > transcripts.fa.transdecoder_dir/longest_orfs.pep.blastp.outfmt6

hmmscan --cpu 8 \
    --domtblout transcripts.fa.transdecoder_dir/longest_orfs.pep.pfam.domtblout \
    $Pfam \
    transcripts.fa.transdecoder_dir/longest_orfs.pep

# （3）进一步筛选ORF
TransDecoder.Predict -t $transcripts \
    --retain_pfam_hits transcripts.fa.transdecoder_dir/longest_orfs.pep.pfam.domtblout \
    --retain_blastp_hits transcripts.fa.transdecoder_dir/longest_orfs.pep.blastp.outfmt6

# ======== Step1 Run blast_hmmer  ========
# blastx -query $transcripts -db TrinotateDB/uniprot_sprot.pep  -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
diamond blastx --query $transcripts --db TrinotateDB/uniprot_sprot.dmnd --threads 24 --max-target-seqs 1 --outfmt 6 > blastx.outfmt6
#blastp -query $proteins -db $uniprot_sprot_db -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
diamond blastp --query $proteins --db $uniprot_sprot_db --threads 24 --max-target-seqs 1 --outfmt 6 > blastp.outfmt6
hmmscan --cpu 32 --domtblout TrinotatePFAM.out $Pfam $proteins > pfam2.log
Trinotate --db TrinotateDB/Trinotate.sqlite --report > myTrinotate.tsv
# 表达矩阵增加注释信息
$HOME/miniforge3/envs/trinity/bin/Trinotate_get_feature_name_encoding_attributes.pl \
    Trinotate_report.xls  > annot_feature_map.txt
${TRINITY_HOME}/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl \
    Trinity_trans.counts.matrix annot_feature_map.txt > Trinity_trans.counts.wAnnot.matrix
Trinotate \
    --db $sqlitedb \
    --init \
    --gene_trans_map Trinity_out/cdhit.fa.gene_trans_map \
    --transcript_fasta transcripts.fa \
    --transdecoder_pep transcripts.fa.transdecoder_dir/transcripts.fa.transdecoder.pep
Trinotate \
    --db $sqlitedb \
    --CPU 32 \
    --gene_trans_map Trinity_out/cdhit.fa.gene_trans_map \
    --transcript_fasta $transcripts \
    --transdecoder_pep $proteins \
    --TrinotateDB TrinotateDB \
    --run "swissprot_blastp swissprot_blastx pfam signalp6 tmhmmv2 infernal EggnogMapper" \
    --use_diamond

Trinotate --db $sqlitedb --LOAD_swissprot_blastp blastp.outfmt6
Trinotate --db $sqlitedb --LOAD_pfam orgdb/TrinotatePFAM.out
# Trinotate --db $sqlitedb --LOAD_signalp <file>
Trinotate --db $sqlitedb --LOAD_EggnogMapper Z100.emapper.annotations
Trinotate --db $sqlitedb --LOAD_swissprot_blastx blastx.outfmt6
# Trinotate --db $sqlitedb --LOAD_infernal <file>

Trinotate \
    --db $sqlitedb \
    --report \
    -E 1e-5 > orgdb/TrinotateReport.xls

download_eggnog_data.py --dbname eggnog --data_dir orgdb
emapper.py \
-i $proteins \
--itype proteins \
--data_dir orgdb \
--output_dir emapper \
--output emapper \
-m diamond \
--cpu 32
