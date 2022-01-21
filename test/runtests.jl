using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using SparseArrays
using Test

#close_chrms = DataFrame(CSV.File("test/close_chroms.csv"))
#far_chrms = DataFrame(CSV.File("test/far_chroms.csv"))


@testset "Simulated Chromosomes"  begin
    include("simpletest.jl")
end


close_chrms = DataFrame(CSV.File("close_chroms.csv"))
far_chrms = DataFrame(CSV.File("far_chroms.csv"))

r = 1500
sig = 750
#r = 3000
#sig = 1500
max_paths = 2

# ldps, allele, gs = DNASeqFISHChromosomeAssignment.find_longest_disjoint_paths(far_chrms, r, sig, max_paths)

res = assign_chromosomes(far_chrms, r, sig, 2, 500)
ldp_final_eq = (res.ldp_allele .== res.final_allele)
ldp_not_classified = (res.ldp_allele .== -1)
allele_not_classified = (res.allele .== -1)

@testset "Far Chromosomes" begin
    @test all(ldp_final_eq .| ldp_not_classified .| allele_not_classified)
    @test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
end

res2 = assign_chromosomes(close_chrms, r, sig, max_paths, 100)
@testset "Close Chromosomes" begin
    @test all(res2.ldp_allele .== res2.final_allele)
    @test all(sort(unique(res2.ldp_allele)) .== [-1, 1, 2])
end
