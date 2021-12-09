using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using SparseArrays
using Test

close_chrms = DataFrame(CSV.File("close_chroms.csv"))
far_chrms = DataFrame(CSV.File("far_chroms.csv"))


r = 1500
sig = 750

dbscan_radius = 1500
min_size = 100
dbscan_min_pnts = 10
max_paths = 2
ldp_overlap_thresh = 0.99



#@testset "Far Chromosomes" begin
res = assign_chromosomes(far_chrms, r, sig, max_paths)
ldp_final_eq = (res.ldp_allele .== res.final_allele)
ldp_not_classified = (res.ldp_allele .== -1)
@test all(ldp_final_eq .| ldp_not_classified)
@test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
#end

#@testset "Close Chromosomes" begin
res = assign_chromosomes(close_chrms, r, sig, max_paths)
@test all(res.ldp_allele .== res.final_allele)
@test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
#end
println("finished testing")
