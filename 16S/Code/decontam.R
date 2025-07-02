# Inspect Library Sizes

# Let’s take a quick first look at the library sizes (i.e. the number of reads) in each sample, as a function of whether that sample was a true positive sample or a negative control:

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

# Identify Contaminants - Frequency
contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
head(contamdf.freq)
table(contamdf.freq$contaminant)
plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="quant_reading") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")