using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using SparseArrays
using Test

#close_chrms = DataFrame(CSV.File("test/close_chroms.csv"))
#far_chrms = DataFrame(CSV.File("test/far_chroms.csv"))
close_chrms = DataFrame(CSV.File("close_chroms.csv"))
far_chrms = DataFrame(CSV.File("far_chroms.csv"))

r = 1500
sig = 750
max_paths = 2

#res = assign_chromosomes(far_chrms, r, sig, max_paths)
res = assign_chromosomes(close_chrms, r, r, sig)
#ldp_final_eq = (res.ldp_allele .== res.final_allele)
dbscan_final_eq = (res.dbscan_allele .== res.final_allele)

#ldp_not_classified = (res.ldp_allele .== -1)

@testset "Far Chromosomes" begin
    #@test all(ldp_final_eq .| ldp_n t_classified)
    #@test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
    @test all(res.dbscan_allele .== res.final_allele)
    @test all(sort(unique(res.final_allele)) .== [-1, 1, 2])
end

#res = assign_chromosomes(close_chrms, r, sig, max_paths)
res = assign_chromosomes(close_chrms, r, r, sig)


@testset "Close Chromosomes" begin
    #@test all(res.ldp_allele .== res.final_allele)
    @test all(res.allele .== res.final_allele)
    #@test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
    @test all(sort(unique(res.final_allele)) .== [-1, 1, 2])
end
