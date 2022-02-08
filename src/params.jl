export ChromSepParams, set_r_dbscan, set_r_ldp, set_r_ldp_nbr, set_min_size, set_dbscan_min_nbrs, set_minProp_unique, set_sigma

"""
    mutable struct ChromSepParams(
        r_dbscan :: Real
        r_ldp :: Real
        r_ldp_nbr :: Real
        sigma :: Real
        min_size :: Int64
        dbscan_min_nbrs :: Int64
        min_prop_unique :: Real
        )
- `r_dbscan` : the DBSCAN search radius
- `r_ldp` : the longest disjoint path search radius
- `r_ldp_nbr` : the search radius to add neighboring loci to LDP clusters after separating a DBSCAN cluster
- `sigma` : defines the weight function for how much to penalize grouping adjacent genomic loci onto the the same chromosome based on their spatial distance from each other.
- `min_size` : the minimum acceptable number of loci on a strand.
- `dbscan_min_pnts` : the minimum allowed number of points for a DBSCAN cluster. DBSCAN clusters with fewer points are discarded.
- `min_prop_unique` : the minimum proportion of unique loci a DBSCAN cluster is allowed to have before considering splitting it into longest disjoint paths.
"""
mutable struct ChromSepParams
    r_dbscan :: Real
    r_ldp :: Real
    r_ldp_nbr :: Real
    sigma :: Real
    min_size :: Int64
    dbscan_min_nbrs :: Int64
    min_prop_unique :: Real
end

"""
    ChromSepParams() 

Constructor for ChromSepParams object using the default parameter values. The user can set the parameter values using parameter setting methods.
"""
function ChromSepParams() 
    return ChromSepParams(500, 700, 500, 250, 30, 12, 0.9)
end

"""
    set_r_dbscan(p :: ChromSepParams, r :: Real)

Set the DBSCAN search radius for clustering loci.
"""
set_r_dbscan(p :: ChromSepParams, r :: Real) = p.r_dbscan = r 

"""
    set_r_ldp(p :: ChromSepParams, r :: Real)

Set the longest disjoint path search radius for clustering loci.
"""
set_r_ldp(p :: ChromSepParams, r:: Real) = p.r_ldp = r 

"""
    set_r_ldp_nbr(p :: ChromSepParams, r :: Real)

Set the search radius for clustering loci near LDPs after separating DBSCAN clusters with LDPs.
"""
set_r_ldp_nbr(p :: ChromSepParams, r :: Real) = p.r_ldp_nbr = r

"""
    set_sigma(p :: ChromSepParams, s :: Real)

Set the sigma for the spatial distance penalty in the longest disjoint paths optimization.
"""
set_sigma(p :: ChromSepParams, s :: Real) = p.sigma = s

"""
    set_min_size(p :: ChromSepParams, min_size :: Real)

Set the minimum size of DBSCAN clusters and minimum length of the longest disjoint paths.
"""
set_min_size(p :: ChromSepParams, min_size :: Real) = p.min_size = min_size

"""
    set_dbscan_min_nbrs(p :: ChromSepParams, mns :: Int64)

Set the minimum number of neighbors a locus must have within the DBSCAN radius to be incuded in a DBSCAN cluster.
"""
set_dbscan_min_nbrs(p :: ChromSepParams, mns :: Int64) = p.dbscan_min_nbrs = mns

"""
    set_minProp_unique(p :: ChromSepParams, mpu :: Real)

If the proportion of unique loci in a DBSCAN cluster is less than this, the algorithm will then check whether the 
average spatial distance between loci of subsequent genomic coordinates is greater than r_ldp. If so, then the DBSCAN cluster
is split into to Longest Disjoint paths.
"""
set_minProp_unique(p :: ChromSepParams, mpu :: Real) = p.min_prop_unique = mpu