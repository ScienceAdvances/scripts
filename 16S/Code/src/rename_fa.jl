
using FASTX,DataFrames,CSV

idmaps = CSV.read("tmp/idmaps.xls", DataFrame, delim="\t")
idmap_dict = Dict(row[1] => row.ASVID for row in eachrow(idmaps))

input = "tmp/dna-sequences.fasta"
output = "Result/02.asv_taxonomy/RepresentSequence.fasta"

reader = open(FASTA.Reader, input)
writer = open(output, "w")
fasta_writer = FASTA.Writer(writer)

for record in reader
    newid = get(idmap_dict, identifier(record), identifier(record))
    write(fasta_writer, FASTA.Record(newid, sequence(record)))
end

close(reader)
close(fasta_writer)
close(writer)
