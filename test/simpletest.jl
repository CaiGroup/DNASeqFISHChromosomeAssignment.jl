using DNASeqFISHChromosomeAssignment
using DataFrames

# make simple test "chromosome"

df = DataFrame(
    Dict(
        "fov" => fill(1,11),
        "cellID" => fill(1,11),
        "chrom" => fill(1,11),
        "x" => fill(5.0,11),
        "y" => Array(10.0:20.0),
        "z" => fill(1.0,11),
        "g" => Array(1:11),
    )
)

res = assign_chromosomes(df, 3, 1.5, 2, 5)