using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using Test

test_files = ["test/E14_r2r1p1c1chr17_r.csv",
"test/E14_r2r1p1c1chr16_r.csv",
#"test/e14r3r2p8c2c2_r.csv",
#"test/e14_r2r1p1c3c7_r.csv",
"test/e14_r2r1p4c6c2_r.csv",
]

test_files = ["E14_r2r1p1c1chr17_r.csv",
"E14_r2r1p1c1chr16_r.csv",
"e14_r2r1p4c6c2_r.csv",
]

#pnts = DataFrame(CSV.File("test/E14_r2r1p1c1chr17.csv"))
#pnts = DataFrame(CSV.File("test/E14_r2r1p1c1chr16.csv"))
#pnts = DataFrame(CSV.File("test/e14r3r2p8c2c2.csv"))
#pnts = DataFrame(CSV.File("test/e14_r2r1p1c3c7.csv"))
#pnts = DataFrame(CSV.File("test/e14_r2r1p4c6c2.csv"))


@testset "Test E14 Data" begin
    for filename in test_files
        println(filename)
        pnts = DataFrame(CSV.File(filename))


        pnts[!,"x"] .*= 103 #nm/pixel
        pnts[!,"y"] .*= 103 #nm/pixel
        pnts[!,"z"] .*= 250 #nm/slice

        ps = pnts[:, ["fov", "cellID", "chrom", "x", "y", "z", "g"]]

        prm = ChromSepParams()

        #res = assign_chromosomes(ps, r_dbscan, r_ldp, s, min_size, r_dbscan, unique_prop_thresh, dbscan_pnts)#, CPLEX.Optimizer)
        res = assign_chromosomes(ps, prm)#, CPLEX.Optimizer)

        #CSV.write(split(filename,",")[1] * "_r.csv", res)
        new_alleles = Array(pnts[:, ["dbscan_allele", "dbscan_ldp_allele", "dbscan_ldp_nbr_allele"]])
        saved_alleles = Array(res[:, ["dbscan_allele", "dbscan_ldp_allele", "dbscan_ldp_nbr_allele"]])
        @test all(new_alleles .== saved_alleles)
    end
end