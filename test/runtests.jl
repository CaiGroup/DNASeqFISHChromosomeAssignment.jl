using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using SparseArrays
using Test
#using CPLEX

#close_chrms = DataFrame(CSV.File("test/close_chroms.csv"))
#far_chrms = DataFrame(CSV.File("test/far_chroms.csv"))
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
#res = assign_chromosomes(far_chrms, r, sig, max_paths)
res = assign_chromosomes(far_chrms, prms)#r, r, sig, r)#, optimizer=CPLEX.Optimizer)
#ldp_final_eq = (res.ldp_allele .== res.final_allele)
dbscan_final_eq = (res.dbscan_allele .== res.dbscan_ldp_nbr_allele)

#ldp_not_classified = (res.ldp_allele .== -1)

@testset "Far Chromosomes" begin
    #@test all(ldp_final_eq .| ldp_n t_classified)
    #@test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
    @test all(res.dbscan_allele .== res.dbscan_ldp_nbr_allele)
    @test all(sort(unique(res.dbscan_ldp_nbr_allele)) .== [-1, 1, 2])
end

#res = assign_chromosomes(close_chrms, r, sig, max_paths)
res = assign_chromosomes(close_chrms, prms)#r, r, sig, r)#, min_size, min_prop_unique, dbscan_min_pnts)#, CPLEX.Optimizer)


@testset "Close Chromosomes" begin
    #@test all(res.ldp_allele .== res.final_allele)
    @test all(res.allele .== res.dbscan_ldp_allele)
    #@test all(sort(unique(res.ldp_allele)) .== [-1, 1, 2])
    @test all(sort(unique(res.dbscan_ldp_nbr_allele)) .== [-1, 1, 2])
end
