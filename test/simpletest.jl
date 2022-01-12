using DNASeqFISHChromosomeAssignment
using DataFrames
using Test

# make simple test "chromosome"

nloci = 20

main = DataFrame(
    Dict(
        "fov" => fill(1,nloci),
        "cellID" => fill(1,nloci),
        "chrom" => fill(1,nloci),
        "x" => fill(5.0,nloci),
        "y" => Array(10.0:(9.0+nloci)),
        "z" => fill(1.0,nloci),
        "g" => Array(1:nloci),
    )
)

res1 = assign_chromosomes(main, 5.0, 1.5, 4, 2, false)

res1_alleles = unique(res1.ldp_allele)
@test length(res1_alleles) == 1 && -1 ∉ res1_alleles

branch = DataFrame(
    Dict(
        "fov" => fill(1,nloci-5),
        "cellID" => fill(1,nloci-5),
        "chrom" => fill(1,nloci-5),
        "x" => Array(6.0:nloci),
        "y" => fill(14, nloci-5),
        "z" => fill(1.0,nloci-5),
        "g" => Array(6:nloci),
    )
)

branched = vcat(main, branch)

res2 = assign_chromosomes(branched, 5.0, 1.5, 4, 2, false)

res2_alleles = unique(res2.ldp_allele)
@test length(res2_alleles) == 1 && -1 ∉ res2_alleles

main2 = copy(main)

main2.x = main2.x .+ 7.5

multi = vcat(branched, main2)

res3 = assign_chromosomes(multi, 5.0, 1.5, 4, 2, false)
