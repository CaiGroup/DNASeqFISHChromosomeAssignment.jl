using DNASeqFISHChromosomeAssignment
using DataFrames
using Test

# make simple test "chromosome"

nloci = 11

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

branched[!,"loc_id"] = fill(1, nrow(branched))
main2[!,"loc_id"] .= fill(2, nrow(main2)) 

multi = vcat(branched, main2)

res3 = assign_chromosomes(multi, 5.0, 1.5, 2, 2, false)

m2_alleles = unique(res3[res3.x .== 12.5, "ldp_allele"])

b_alleles = unique(res3[res3.x .!= 12.5, "ldp_allele"])

@test length(m2_alleles) == 1 && -1 ∉ m2_alleles
@test length(b_alleles) == 1 && -1 ∉ b_alleles
@test m2_alleles[1] != b_alleles[1]

sc2 = copy(branch)

sc2[!,"loc_id"] = fill(3, nrow(sc2))

sc2.y .-= 50

multi2 = vcat(multi,sc2)

res4 = assign_chromosomes(multi2, 5.0, 1.5, 2, 2, false)

@test all(res4.loc_id .== res4.ldp_allele)
