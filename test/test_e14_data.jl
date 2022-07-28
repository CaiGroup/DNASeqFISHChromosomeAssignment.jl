using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using Test
using GLPK


test_files = [
"E14_r2r1p1c1chr16.csv", # far apart chromsomes
"e14_r2r1p4c6c2.csv" #close together chromosomes
]


@testset "Test E14 Data" begin
    for filename in test_files
        pnts = DataFrame(CSV.File(filename))

        #pnts[!,"x"] .*= 103 #nm/pixel
        #pnts[!,"y"] .*= 103 #nm/pixel
        #pnts[!,"z"] .*= 250 #nm/slice

        ps = pnts[:, ["fov", "cellID", "chrom", "x", "y", "z", "g"]]

        prm = ChromSepParams()
        set_dbscan_min_nbrs(prm, 6)

        res = assign_chromosomes(ps, prm, GLPK.Optimizer, false)

        new_alleles = Array(pnts[:, ["dbscan_allele", "dbscan_ldp_allele", "dbscan_ldp_nbr_allele"]])
        saved_alleles = Array(res[:, ["dbscan_allele", "dbscan_ldp_allele", "dbscan_ldp_nbr_allele"]])
        @test all(new_alleles .== saved_alleles)
    end
end