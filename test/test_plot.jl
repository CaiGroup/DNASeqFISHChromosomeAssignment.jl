using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment

df = DataFrame(CSV.File("ab4.csv"))
df[!,"children"] = eval.(df[!,"children"])
df[!,"parents"] = eval.(df[!,"parents"])


DNASeqFISHChromosomeAssignment.plot_chrom(df, 28, 5)