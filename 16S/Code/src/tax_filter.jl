using FASTX,DataFrames,CSV,XLSX

rdp = CSV.read("Result/02.asv_taxonomy/rdp_taxonomy.xls", DataFrame, delim="\t",header=false)
silva = CSV.read("Result/02.asv_taxonomy/silva_taxonomy.xls", DataFrame, delim="\t")

asv = intersect(rdp.Column1, silva.ASVID)

keep_list = []
for i in asv
    x=rdp[rdp.Column1.==i,:Column4]
    x2=split(x[1], ",")
    x3=[x[2] for x in split.(x2,":")]
    x4=replace.(x3,"/" => "-")
    x_pad = vcat(x4, fill("", 7 - length(x4)))

    y=silva[silva.ASVID.==i,2:8]
    y2=collect(y[1, :])
    y3 = [ismissing(val) ? "" : val for val in y2]
    y4 = [length(x) >= 2 ? x[2] : "" for x in split.(y3, "__")]
    y_pad = vcat(y4, fill("", 7 - length(y4)))
    
    if sum(x_pad .== y_pad) > 2
        push!(keep_list, i)
    end
end

df=silva[[x ∈ keep_list for x in silva.ASVID] ,:];
XLSX.writetable("Result/02.asv_taxonomy/taxonomy.xlsx", df)
CSV.write("Result/02.asv_taxonomy/taxonomy.xls",df, delim="\t")

ASV = CSV.read("Result/02.asv_taxonomy/ASV.xls", DataFrame, delim="\t")
df=ASV[[x ∈ keep_list for x in ASV.ASVID] ,:];
XLSX.writetable("Result/02.asv_taxonomy/ASV.xlsx", df)
rename!(df, names(df)[1] => "#OTU");
CSV.write("Result/02.asv_taxonomy/ASV.xls",df, delim="\t");
