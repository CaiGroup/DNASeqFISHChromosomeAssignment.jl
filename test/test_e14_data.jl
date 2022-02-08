using CSV
using DataFrames
using DNASeqFISHChromosomeAssignment
using Test

test_files = ["E14_r2r1p1c1chr17_r.csv",
"E14_r2r1p1c1chr16_r.csv",
"e14_r2r1p4c6c2_r.csv",
]


@testset "Test E14 Data" begin
    for filename in test_files
        pnts = DataFrame(CSV.File(filename))

        pnts[!,"x"] .*= 103 #nm/pixel
        pnts[!,"y"] .*= 103 #nm/pixel
        pnts[!,"z"] .*= 250 #nm/slice

        ps = pnts[:, ["fov", "cellID", "chrom", "x", "y", "z", "g"]]

        prm = ChromSepParams()

        res = assign_chromosomes(ps, prm)

        new_alleles = Array(pnts[:, ["dbscan_allele", "dbscan_ldp_allele", "dbscan_ldp_nbr_allele"]])
        saved_alleles = Array(res[:, ["dbscan_allele", "dbscan_ldp_allele", "dbscan_ldp_nbr_allele"]])
        @test all(new_alleles .== saved_alleles)
    end
end