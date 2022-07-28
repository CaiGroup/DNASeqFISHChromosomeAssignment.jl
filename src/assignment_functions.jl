using DataFrames
using JuMP
using NearestNeighbors
using Graphs
using SparseArrays
using GLPK
using Clustering
using Statistics

"""
	assign_chromosomes(pnts :: DataFrame, params :: ChromSepParams, optimizer)

This function assigns genomic loci from DNA seqFISH+ data to indiviual chromosomes using DBSCAN and Longest Disjoint Paths (LDP). First, genomic loci
are clustering using DBSCAN. Then DBSCAN clusters that contain a proportion of unique genomic loci less than the min_prop_unique parameter set in the params
object and whose average spatial distance between subsequent spatial loci is greater than the r_ldp parameter (also et in the params object) will be split,
using the LDP algorithm. If these conditions are not met, then we find one LDP that is a subset of the cluster.

LDP works by finding the longest one or two weighted paths of loci of increasing genomic coordinate within the a DBSCAN cluster. More specifically, it builds a directed acyclic graph where every locus has an out going edge to genomically
down stream loci within the search radius. Edges are weighted by:

``
w = \\frac{1}{g_p - g_c} × e^{-\\frac{||x_p - x_c||^2}{2σ^2}}
``

where x denotes a vector of spatial coordinates, g denotes genomic coordinates, and p and c subscripts denote parent and child nodes respectively.
LDP then uses Integer Programming to find the one or two disjoint paths that maximize the sum of
the weights of all edges in the paths.

Arguments
- df: Takes a DataFrame (pnts) of decoded DNA seqFISH data where each row describes a locus. Must include columns:
	- fov : the field of view in which locus was imaged
	- cellID : the cell in which the locus was imaged
	- chrom : chromosome on which locus is located
	- x : x spatial coordinate of locus
	- y : y spatial coordinate of locus
	- z : z spatial coordinate of locus
	- g : genomic coordinate of locus
- prms : ChromSepParams object containing parameters.
- optimizer : optional - JuMP compatible optimizer for performing Longest Disjoint Paths optimization. Uses GLPK by default.


returns a copy of the input points sorted by fov, cellID, chrom, and g with a
new columns: 
	`dbscan_allele` : DBSCAN clusters of each locus
	`dbscan_ldp_allele` : ldp subclusters of each dbscan cluster
	`dbscan_ldp_nbr_allele` : dbscan clusters that are not split, and expands split ldp clusters by assigning ldp neighbors to their cluster.


allele = -1 means that a locus was not assigned to a chromosome.
allele > 0 indicates the allele number of the chromosome that the locus was assigned to.

"""
function assign_chromosomes(pnts :: DataFrame,
	 						prms :: ChromSepParams,
							optimizer = GLPK.Optimizer,
							auto_set_rs = true)

	@assert prms.min_size > 1
	sort!(pnts, [:fov, :cellID, :chrom, :g])
	chrms = groupby(pnts,[:fov, :cellID, :chrom])

	# first assign chromosomes by DBSCAN
	assn_chrms_dbscan(chrm) = assign_loci_dbscan(chrm, prms, auto_set_rs)
	dbscan_res = transform(assn_chrms_dbscan, chrms)
	println("finished DBSCAN")

	dbscan_grps = groupby(dbscan_res, [:fov, :cellID, :chrom, :dbscan_allele])
	
	if auto_set_rs
		unambiguous_chrms = filter(chrm -> chrm.dbscan_allele[1] != -1 && length(unique(chrm.g))/nrow(chrm) > prms.min_prop_unique, dbscan_grps)

		dists=vcat(map(get_nbr_loci_dists, collect(unambiguous_chrms))...)

		#unambiguous_res = transform(assn_chrms, unambiguous_chrms)
		r_ldp = mean(dists)+std(dists)
		sigma = r_ldp/2
		set_r_ldp(prms, r_ldp)
		set_sigma(prms, sigma)
	end

	assn_chrms_preset_r(chrm) = assign_loci_ldp(chrm, prms, optimizer)

	res = transform(assn_chrms_preset_r, dbscan_grps)
	return res
end

function assign_loci_dbscan(chrm, prms :: ChromSepParams, auto_choose_r=true)
	if nrow(chrm) < prms.min_size || nrow(chrm) <= prms.dbscan_min_nbrs+1 || nrow(chrm) < 4
		return DataFrame(Dict("dbscan_allele"=>fill(-1, nrow(chrm))))
	end
	if auto_choose_r
		set_r_dbscan(prms, auto_choose_dbscan_r(chrm, prms))
	end
	dbscan_clusters = cluster_chromosomes_DBSCAN(chrm, prms)
	dbscan_allele = get_allele_col(chrm, dbscan_clusters)
	return DataFrame(Dict("dbscan_allele"=>dbscan_allele))
end

function assign_loci_ldp(chrm, prms :: ChromSepParams, optimizer=GLPK.Optimizer)
	#if auto_choose_r
	#	set_r_dbscan(prms, auto_choose_dbscan_r(chrm, prms))
	#end
	#dbscan_clusters = cluster_chromosomes_DBSCAN(chrm, prms)
	#dbscan_allele = get_allele_col(chrm, dbscan_clusters)
	#dbscan_ldp_allele, dbscan_ldp_nbr_allele = get_DBSCAN_cluster_LDPs(chrm, dbscan_clusters, dbscan_allele, prms, optimizer)
	if chrm.dbscan_allele[1] == -1
		return DataFrame(Dict("dbscan_ldp_allele"=>fill(-1, nrow(chrm)), "dbscan_ldp_nbr_allele"=>fill(-1, nrow(chrm))))
	else
		dbscan_ldp_allele, dbscan_ldp_nbr_allele = get_DBSCAN_cluster_LDPs(chrm, prms, optimizer)
		return_df = DataFrame(Dict("dbscan_ldp_allele"=>dbscan_ldp_allele, "dbscan_ldp_nbr_allele"=>dbscan_ldp_nbr_allele))
		return return_df
	end
end

function auto_choose_dbscan_r(chrm, prms :: ChromSepParams)
	coords = Array([chrm.x chrm.y chrm.z]')
    tree = KDTree(coords)
	idxs, dists = knn(tree, coords, prms.dbscan_min_nbrs+1)
	min_pnt_core_rads = maximum.(dists)
	r_dbscan = mean(min_pnt_core_rads) + std(min_pnt_core_rads)
	#println("mean(dists) ", mean(min_pnt_core_rads), " sd(dists) ", std(min_pnt_core_rads))
	#println("r_dbscan $r_dbscan")
	return r_dbscan
end

function auto_choose_ldp_r_sigma(chrm, params)
	g_grpd = groupby(chrm,:g)
	#g_grpd_a = Array(g_grpd)
	locus_copynum = combine(nrow, g_grpd)
	# if below uniqueness threshold (two chromosomes)
	if nrow(locus_copynum)/nrow(chrm) < params.min_prop_unique
		println("min_prop_unique: ", params.min_prop_unique)
		println("prop unique: ", nrow(locus_copynum)/nrow(chrm))
		#copynum_dict = Dict([(row.g, row.nrow) for row in eachrow(locus_copynum)])
		dists = Float64[]
		sizehint!(dists, 2*sum(locus_copynum.nrow .> 1))
		for (locus_ind, locus_df) in enumerate(g_grpd[1:end-1])
			if nrow(g_grpd[locus_ind+1]) > 1 
				for locus in eachrow(DataFrame(locus_df))
					push!(dists, get_min_nbr_dist(g_grpd[locus_ind+1], locus))
				end
			end
			if nrow(locus_df) > 1
				for nbr in eachrow(DataFrame(g_grpd[locus_ind+1]))
					push!(dists, get_min_nbr_dist(locus_df, nbr))
				end
			end
		end 
	else
		filtered_chrm = DataFrame(filter(grp -> nrow(grp) == 1, g_grpd))
		#println(filtered_chrm)
		nrows = nrow(filtered_chrm)
    	dists = (Array(filtered_chrm[2:end,["x","y","z"]] .- filtered_chrm[1:(nrows-1),["x","y","z"]])).^2
		dists = sqrt.(sum(dists, dims=2))
	end
	println("mean(dists) ", mean(dists), " std(dists) ", std(dists))
	r_ldp = mean(dists) + std(dists)
	sigma = r_ldp/2
	println("r_ldp $r_ldp, sigma $sigma")
	return r_ldp, sigma
end

# find rs from dbscan separated chromsomes
function find_reasonable_rs(unambiguous_chrms)
	groupby(unambiguous_chrms, [:])
	
end

function get_min_nbr_dist(df1, df2)
	sqd = (DataFrame(df1)[:, ["x","y","z"]] .- DataFrame(df2)[:, ["x","y","z"]]).^2
	return minimum(sqrt.(sum(Array(sqd),dims=2)))
end

function find_longest_disjoint_paths(chrm, prms :: ChromSepParams, max_strands :: Int64, optimizer)
	A, W, g = get_neighbors(chrm, prms.r_ldp, prms.sigma)
	R, WR = get_trans_red_w(chrm, prms.r_ldp, prms.sigma, A, W, g)
	g2, W2 = rule_out_edges(A, W, WR)
	return optimize_paths(chrm, g2, W2, prms.min_size, max_strands, optimizer)
end

function get_neighbors(chr, r, sig)
    nloci = size(chr)[1]
    chr = DataFrame(chr)
    sort!(chr, "g")

    # get spatial neighbors
    coords = Array([chr.x chr.y chr.z]')
    tree = KDTree(coords)
    nbrs = inrange(tree, coords, r)

    # build neighbor adjacency graph weighted by genomic distance
    A = spzeros(Bool, nloci, nloci)
    W = spzeros(Float64, nloci, nloci)
    g = SimpleDiGraph(nloci)

    for (locus, nbrs) in enumerate(nbrs)
        # add outgoing edges to genomically down stream spatial neighbors
        for nbr in nbrs
            if chr[nbr,"g"] > chr[locus,"g"]
				genomic_weight = 1/(abs(chr[nbr,"g"] - chr[locus,"g"]))
                spatial_weight = exp(-((chr[locus,"x"] - chr[nbr,"x"])^2 +
                                       (chr[locus,"y"] - chr[nbr,"y"])^2 +
                                       (chr[locus,"z"] - chr[nbr,"z"])^2)/(2 * sig^2))
                A[locus, nbr] = 1
                W[locus, nbr] = genomic_weight * spatial_weight
                add_edge!(g, locus, nbr)
            end
        end
    end

    return A, W, g
end

function get_trans_red_w(chr, r, sig, A, W, g)
    # compute transitive reduction
    B = adjacency_matrix(transitiveclosure(g))
    R = A .* (0 .== (A*B))

    #get weights of transitive reduction
    WR = R .* W
    return R, WR
end

function rule_out_edges(A, W, WR)
    nnodes = size(A)[1]
	A2 = spzeros(nnodes, nnodes)
    W2 = spzeros(size(A)...)
    for i in 1:nnodes
        min_weight_reduced_arc = isempty(nonzeros(WR[i,:])) ? 0.0 : minimum(nonzeros(WR[i,:]))
        outgoing_edges_to_keep = (W[i,:] .>= min_weight_reduced_arc)
        W2[i,outgoing_edges_to_keep] .= W[i,outgoing_edges_to_keep]
    end
    A2[1:nnodes, 1:nnodes] = W2 .> 0
    g2 = DiGraph(A2)
    return g2, W2
end

function optimize_paths(chrm, g :: DiGraph, W :: SparseMatrixCSC, min_size :: Int64, max_strands :: Int64, optimizer)
	nloci = nv(g)
	src = nloci+1
	dst = src+1
	locus_edges =Tuple.(collect(Graphs.edges(g)))
	add_vertices!(g, 2)
	for i in 1:nloci
		add_edge!(g, src, i)
		add_edge!(g, i, dst)
	end
	model = Model(optimizer)

	g_edges = Tuple.(collect(Graphs.edges(g)))

	@variable(model, x[g_edges], Bin, container=SparseAxisArray)
	@constraint(model, sum(x[(src,nbr)] for nbr in 1:nloci) <= max_strands)
	@constraint(model, sum(x[(nbr,dst)] for nbr in 1:nloci) <= max_strands)
	@constraint(model, [i = 1:nloci], sum(x[(nbr,i)] for nbr in inneighbors(g, i)) == sum(x[(i,nbr)] for nbr in outneighbors(g, i)))
	@constraint(model, [i = 1:nloci], sum(x[(nbr,i)] for nbr in inneighbors(g, i)) <= 1)

	@objective(model, Max, sum(x[e]*W[e...] for e in locus_edges))

	optimize!(model)

	evals = value.(x)
	gres = DiGraph(nloci)
	for e in locus_edges
		evals[e] == 1 ? add_edge!(gres, e[1], e[2]) : nothing
	end

	wccs = weakly_connected_components(gres)

	lccs = length.(wccs)
	ldps = wccs[lccs .>= min_size]

	allele = get_allele_col(chrm, ldps)

	return ldps, allele
end

function get_allele_col(chrm, grps)
	allele = fill(-1, nrow(chrm))
	for (i, grp) in enumerate(grps)
		allele[grp] .= i
	end
	return allele
end

function cluster_chromosomes_DBSCAN(chrms, prms)
	if nrow(chrms) <= 3
		return []
	end
	points = Array([chrms.x chrms.y chrms.z]')
	dbr = dbscan(points, prms.r_dbscan, min_neighbors=prms.dbscan_min_nbrs, min_cluster_size=prms.min_size)
	dbscan_clusters = [sort(vcat(dbc.core_indices,  dbc.boundary_indices)) for dbc in dbr]
	return dbscan_clusters
end


"""
Get a strict Longest dijsoint path or two for every DBSCAN cluster. When two LDPs are found, also group points near each ldp

"""
#ToDo: Need to figure out how to deal with subsetting multiple dbscan clusters of same chromosome

#function get_DBSCAN_cluster_LDPs(chrm, dbscan_clusters, dbscan_allele, prm, optimizer)
function get_DBSCAN_cluster_LDPs(chrm, prm, optimizer)
	dbscan_ldp_allele = fill(-1, nrow(chrm))
	dbscan_ldp_nbr_allele = fill(-1, nrow(chrm))
	n_new_alleles = 0
	ldp_allele_nums = [-1, chrm.dbscan_allele[1]] #vcat([-1], Array(1:length(dbscan_clusters)))

	#for (i, c) in enumerate(dbscan_clusters)
		#if auto_set_r
			#r_ldp, sigma = auto_choose_ldp_r_sigma(chrm[c,:], prm)
			#set_r_ldp(prm, r_ldp)
			#set_sigma(prm, sigma)
		#end
	if length(unique(chrm.g))/length(chrm.g) < prm.min_prop_unique && mean_g_nbr_spat_dist(chrm) > prm.r_ldp
		ldps, dbscan_ldp_allele = find_longest_disjoint_paths(chrm, prm, 2, optimizer)
		if length(ldps) > 0
			dbscan_ldp_allele[ldps[1]] .= 1
			dbscan_ldp_nbr_allele[ldps[1]] .= 1
			if length(ldps) == 2
				n_new_alleles += 1
				#new_allele= length(dbscan_clusters) + n_new_alleles
				new_allele= 1 + n_new_alleles
				push!(ldp_allele_nums, new_allele)
				dbscan_ldp_allele[ldps[2]] .= new_allele
				dbscan_ldp_nbr_allele[ldps[2]] .= new_allele
			else
				dbscan_ldp_nbr_allele .= chrm.dbscan_allele
			end
		end
	else
		dbscan_ldp_nbr_allele .= chrm.dbscan_allele
		ldp, dbscan_c_ldp_allele = find_longest_disjoint_paths(chrm, prm, 1, optimizer)
		if length(ldp) > 0
			dbscan_ldp_allele[ldp[1]] .= chrm.dbscan_allele[1]
		end
	end
	dbscan_ldp_nbr_allele = assign_ldp_neighbors(chrm, dbscan_ldp_nbr_allele, prm.r_ldp_nbr, ldp_allele_nums)
	return dbscan_ldp_allele, dbscan_ldp_nbr_allele
end

function get_nbr_loci_dists(loci)
	println("nrow(loci): ", nrow(loci))
	println("dbscan_allele: ", loci.dbscan_allele[1])
	if nrow(loci) < 2
		return []
	end
	dists = Array(loci[2:end,["x","y","z"]] .- loci[1:end-1,["x","y","z"]])
	dists .^= 2
	return sqrt.(sum(dists, dims=2))
end

function mean_g_nbr_spat_dist(loci)
	#nrows = nrow(loci)
    #dists = Array(loci[2:end,["x","y","z"]] .- loci[1:(nrows-1),["x","y","z"]])
	#dists .^= 2
	#mean_dist = mean(sqrt.(sum(dists, dims=2)))
	return mean(get_nbr_loci_dists(loci)) #mean_dist
end

"""
The longest disjoint leaves many points unassigned that do not fit into the path. These points may be replication sites, or
otherwise legitimate. This function counts how many points within a given radius of each unassigned locus are in each
of the longest disjoint paths, and assigned the points to the same allele as the longest disjoint path with the most
loci in that radius.
"""
function assign_ldp_neighbors(pnts, dbscan_ldp_nbr_allele, r, ldp_allele_nums)
	# make KDTrees of loci assigned to each allele
	cluster_nums = sort(filter(c -> c != -1, unique(dbscan_ldp_nbr_allele)))
	c_trees = [KDTree(Array(Array(pnts[dbscan_ldp_nbr_allele .== c,["x","y","z"]])')) for c in cluster_nums]
	unassigned_loci = filter(locus -> dbscan_ldp_nbr_allele[locus] == -1, Array(1:nrow(pnts)))
	unassigned_coords = Array(Array(pnts[dbscan_ldp_nbr_allele .== -1, ["x","y","z"]])')
	n_allele_nbrs = zeros(sum(dbscan_ldp_nbr_allele .== -1), length(ldp_allele_nums))
	for (i, c_tree) in enumerate(c_trees)
		n_allele_nbrs[:, i + 1] .= length.(inrange(c_tree, unassigned_coords, r))
	end

	ldp_nbr_allele = broadcast(x -> ldp_allele_nums[x], getindex.(argmax(n_allele_nbrs, dims=2),2))
	ldp_nbr_allele = reshape(ldp_nbr_allele, length(ldp_nbr_allele))
	
	dbscan_ldp_nbr_allele[unassigned_loci] .= ldp_nbr_allele 
	return dbscan_ldp_nbr_allele
end