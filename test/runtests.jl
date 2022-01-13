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

#r = 1500
#sig = 750
r = 3000
sig = 1500
max_paths = 2

res = assign_chromosomes(far_chrms, r, sig, max_paths)
ldp_final_eq = (res.ldp_allele .== res.final_allele)
ldp_not_classified = (res.ldp_allele .== -1)

@testset "Far Chromosomes" begin
    @test all(ldp_final_eq .| ldp_not_classified)
    @test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
end

res = assign_chromosomes(close_chrms, r, sig, max_paths)
@testset "Close Chromosomes" begin
    @test all(res.ldp_allele .== res.final_allele)
    @test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
end
