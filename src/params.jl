export ChromSepParams, set_r_dbscan, set_r_ldp, set_r_ldp_nbr, set_min_size, set_dbscan_min_pnts, set_minProp_unique, set_sigma


mutable struct ChromSepParams
    r_dbscan :: Real
    r_ldp :: Real
    r_ldp_nbr :: Real
    sigma :: Real
    min_size :: Int64
    dbscan_min_pnts :: Int64
    min_prop_unique :: Real
end

function ChromSepParams() 
    return ChromSepParams(500, 700, 500, 250, 30, 12, 0.9)
end

set_r_dbscan(p :: ChromSepParams, r :: Real) = p.r_dbscan = r 

set_r_ldp(p :: ChromSepParams, r:: Real) = p.r_ldp = r 

set_r_ldp_nbr(p :: ChromSepParams, r :: Real) = p.r_ldp_nbr = r

set_sigma(p :: ChromSepParams, s :: Real) = p.sigma = s

set_min_size(p :: ChromSepParams, min_size :: Real) = p.min_size = min_size

set_dbscan_min_pnts(p :: ChromSepParams, mps :: Int64) = p.dbscan_min_pnts = mps

set_r_ldp_nbr(p :: ChromSepParams, s :: Real) = p.sigma = s


set_minProp_unique(p :: ChromSepParams, mpu :: Real) = p.min_prop_unique = mpu