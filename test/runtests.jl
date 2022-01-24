using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using SparseArrays
using Test

close_chrms = DataFrame(CSV.File("close_chroms.csv"))
far_chrms = DataFrame(CSV.File("far_chroms.csv"))

dbscan_r_min = 400
dbscan_r_max = 800
dbscan_r_inc = 50
r_ldp = 800
sig = 300
overlap_thresh = 0.1

#res = assign_chromosomes(far_chrms, r, sig, max_paths)
res = assign_chromosomes(far_chrms, dbscan_r_min, dbscan_r_max, dbscan_r_inc, r_ldp, sig, overlap_thresh)
#ldp_final_eq = (res.ldp_allele .== res.final_allele)
ldp_not_classified = (res.allele .== -1)

@testset "Far Chromosomes" begin
    #@test all(ldp_final_eq .| ldp_not_classified)
    @test all((res.allele .== res.final_allele) .| (res.allele .== -1) .| (res.final_allele .== -1))
    @test all(sort(unique(res.final_allele)) .== [-1, 1, 2])
end

c30chr15 = DataFrame(CSV.File("NMuMGr1_c30c15.csv"))
res4 = assign_chromosomes(c30chr15, r, r, sig,30, 0.9)
all(c30chr7[:, "final_allele"] .== c30chr7[:, "dbscan_allele"])


c30chr7 = DataFrame(CSV.File("NMuMGr1_c30chr7.csv"))
res3 = assign_chromosomes(c30chr7, r, r, sig, max_paths,30, 0.9)
all(c30chr7[:, "final_allele"] .== c30chr7[:, "dbscan_allele"])

res = assign_chromosomes(close_chrms, r, sig, max_paths)
@testset "Close Chromosomes" begin
    #@test all(res.ldp_allele .== res.final_allele)
    @test all(res.allele .== res.final_allele .| res.allele .== -1 .| res.final_allele .== -1)
    @test all(sort(unique(res.final_allele)) .== [-1, 1, 2])
end

