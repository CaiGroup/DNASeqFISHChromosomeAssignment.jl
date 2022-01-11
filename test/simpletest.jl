using DNASeqFISHChromosomeAssignment
using DataFrames
using Test

# make simple test "chromosome"

main = DataFrame(
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

res1 = assign_chromosomes(main, 5.0, 1.5, 2, 2, false)

res1_alleles = unique(res1.ldp_allele)
@test length(res1_alleles) == 1 && -1 ∉ res1_alleles

branch = DataFrame(
    Dict(
        "fov" => fill(1,6),
        "cellID" => fill(1,6),
        "chrom" => fill(1,6),
        "x" => Array(6.0:11.0),
        "y" => fill(14, 6),
        "z" => fill(1.0,6),
        "g" => Array(6:11),
    )
)

branched = vcat(main, branch)

res2 = assign_chromosomes(branched, 5.0, 1.5, 2, 2, false)
