using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using SparseArrays
using Test

close_chrms = DataFrame(CSV.File("close_chroms.csv"))
far_chrms = DataFrame(CSV.File("far_chroms.csv"))

r = 1500
sig = 750
max_paths = 2
min_size=100
min_prop_unique=0.9
dbscan_min_pnts = 100
prms = ChromSepParams()
set_min_size(prms, 100)
set_r_ldp(prms, 1500)
set_r_dbscan(prms,1500)
set_sigma(prms, 750)

res = assign_chromosomes(far_chrms, prms)

@testset "Far Chromosomes" begin
    @test all(res.dbscan_allele .== res.dbscan_ldp_nbr_allele)
    @test all(sort(unique(res.dbscan_ldp_nbr_allele)) .== [-1, 1, 2])
end

res = assign_chromosomes(close_chrms[:, ["fov", "cellID", "chrom", "x", "y", "z", "g"]], prms)

@testset "Close Chromosomes" begin
    @test all(close_chrms.dbscan_ldp_allele .== res.dbscan_ldp_allele)
    @test all(sort(unique(res.dbscan_ldp_nbr_allele)) .== [-1, 1, 2])
end

include("test_e14_data.jl")