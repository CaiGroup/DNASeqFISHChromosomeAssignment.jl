using DataFrames
using JuMP
using NearestNeighbors
using Graphs
using SparseArrays
using GLPK
using Clustering

"""
	assign_chromosomes(pnts :: DataFrame,
	 			search_radius :: Real,
				sigma :: Real,
				max_strands :: Int64,
				min_size :: Int64 = 2,
				dbscan :: Bool = true,
				dbscan_min_pnts :: Int64 = 10,
				ldp_overlap_thresh :: Float64 = 0.99)

This function assigns genomic loci from DNA seqFISH+ data to indiviual chromosomes using DBSCAN and Longest Disjoint Paths (LDP).
LDP includes spatial and genomic information to effectely separate chromosome copies that are close to each other
in the cell, while DBSCAN is able to classify repiclated loci onto chromosomes.


LDP builds a directed acyclic graph where every locus has an out going edge to genomically
down stream loci within the search radius. Edges are weighted by:

``
w = (g_{max} - (g_p - g_c)) × e^{-\\frac{||x_p - x_c||^2}{2σ^2}}
``

where x denotes a vector of spatial coordinates, g denotes genomic coordinates, and p and c subscripts denote parent and child nodes respectively. ``g_{max}``
is the maximum genomic coordinate for that chromosome found in the cell.

LDP finds up to max_strands disjoint paths that maximize the sum of
the weights of all edges in the paths.

Arguments
- df: Takes a DataFrame (pnts) of decoded DNA seqFISH data where each row describes a locus. Must include columns:
	- fov : the field of view in which locus was imaged
	- cellID : the cell in which the locus was imaged
	- chrom : chromosome on which locus is located
	- x : x spatial coordinate of locus (Must be Floats)
	- y : y spatial coordinate of locus (Must be Floats)
	- z : z spatial coordinate of locus (Must be Floats)
	- g : genomic coordinate of locus
- `search_radius` : gives how far in space to search for neighboring loci on the same
chromosome strand for each locus.
- `sigma` : defines the weight function for how much to penalize grouping adjacent genomic loci onto the the same chromosome based on their spatial distance from each other.
- `max_strands` : the maximum number of copies of each chromosome that we expect to find.
- `min_size` : the minimum acceptable number of loci on a strand.
- `dbscan` : if true, run both DBSCAN and LDP, then find the consensus, else use only LDP.
- `dbscan_min_pnts` : the minimum allowed number of points for a DBSCAN cluster. DBSCAN clusters with fewer points are discarded .
- `ldp_overlap_thresh` : the minimum proportion of LDP assigned loci that DBSCAN assigned the same for the two clusters two be merged

returns a copy of the input points sorted by fov, cellID, chrom, and g with a
new column: `ldp_allele`. If the `dbscan` argument is true, also includes `dbscan_allele` and `final_allele` columns.

allele = -1 means that a locus was not assigned to a chromosome.
allele > 0 indicates the chromsome number that the locus was assigned to.

If the size of an intersection between a DBSCAN and LDP cluster divided by the size of the LDP cluster is more than than `ldp_overlap_thresh`, we say that the two clusters are the same.

If every LDP cluster overlaps with exactly one DBSCAN cluster, and Every DBSCAN cluster overlaps with one or zero LDP clusters, we say the methods agreed, and the final clusters are the
unions of the overlapping LDP and DBSCAN clusters. Otherwise, the final clusters and the LDP clusters.

In the following example, LDP detects that two copies of a chromosome in the `close_chroms.csv` test data are distinct even though they are spatially close.

```
julia> first(data, 5)
5×7 DataFrame
 Row │ fov    cellID  chrom    x        y          z        g
     │ Int64  Int64   String7  Float64  Float64    Float64  Int64
─────┼────────────────────────────────────────────────────────────
   1 │     0      28  chr6     94128.3  1.23568e5  1631.25      1
   2 │     0      28  chr6     92462.6  1.23556e5  2808.75      6
   3 │     0      28  chr6     94006.9  1.23476e5  1833.25      8
   4 │     0      28  chr6     93877.1  1.23784e5  2195.5      10
   5 │     0      28  chr6     92537.5  1.23505e5  2940.0      11

julia> r = 1500
julia> sig = 750
julia> max_paths = 2

julia> res = assign_chromosomes(data, r, sig, max_paths)

julia> first(res,5)
5×10 DataFrame
 Row │ fov    cellID  chrom    x        y          z        g    ldp_allele  dbscan_allele  final_allele
     │ Int64  Int64   String7  Float64  Float64    Float64  Int64  Int64       Int64          Int64
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │     0      28  chr6     94128.3  1.23568e5  1631.25      1           1              1             1
   2 │     0      28  chr6     92462.6  1.23556e5  2808.75      6           2              1             2
   3 │     0      28  chr6     94006.9  1.23476e5  1833.25      8           1              1             1
   4 │     0      28  chr6     93877.1  1.23784e5  2195.5      10           1              1             1
   5 │     0      28  chr6     92537.5  1.23505e5  2940.0      11           2              1             2

```
"""
function assign_chromosomes(pnts :: DataFrame,
	 						search_radius :: Real,
							sigma :: Real,
							max_strands :: Int64,
							min_size :: Int64 = 100,
							dbscan :: Bool = true,
							dbscan_min_pnts :: Int64 = 10,
							ldp_overlap_thresh :: Float64 = 0.99)

	@assert min_size > 1
	if nrow(pnts) > 1
		@assert typeof(pnts.x) == Vector{Float64}
		@assert typeof(pnts.y) == Vector{Float64}
		@assert typeof(pnts.z) == Vector{Float64}
	end

	sort!(pnts, [:fov, :cellID, :chrom, :g])
	chrms = groupby(pnts,[:fov, :cellID, :chrom])
	assn_chrms(chrm) = assign_loci(chrm, search_radius, sigma, max_strands, min_size, dbscan, dbscan_min_pnts, ldp_overlap_thresh)
	res = transform(assn_chrms, chrms)
	#renamed_col = combine(res, Not(:x1), :x1 => :allele)
	return res #renamed_col
end

function assign_loci(chrm, r :: Real, sig :: Real, max_strands :: Int64, min_size :: Int64, dbscan, dbscan_min_pnts, ldp_overlap_thresh)
	#chrm = sort(_chrm, :g)
	println("fov: ", chrm[1, "fov"], ", cell: ", chrm[1, "cellID"], ", ", chrm[1,"chrom"])
	ldps, ldp_allele = find_longest_disjoint_paths(chrm, r, sig, max_strands, min_size)
	return_df = DataFrame(Dict("ldp_allele"=>ldp_allele))
	if dbscan
		dbscan_clusters = cluster_chromosomes_DBSCAN(chrm, r, dbscan_min_pnts, min_size)
		return_df[!, "dbscan_allele"] = get_allele_col(chrm, dbscan_clusters)
		return_df[!, "final_allele"] = compare_LDP_DBSCAN(ldps, dbscan_clusters, chrm, max_strands,ldp_overlap_thresh)
	end
	return return_df
end

function find_longest_disjoint_paths(chrm, r :: Real, sig :: Real, max_strands :: Int64, min_size :: Int64)
	A, W, g = get_neighbors(chrm, r, sig)
	R, WR = get_trans_red_w(chrm, r, sig, A, W, g)
	g2, W2 = rule_out_edges(A, W, WR)
	#wcc_gs, wcc_Ws, wccs = get_connected_components(g2, W2)
	wccs = weakly_connected_components(g2)
	allele = fill(-1, nrow(chrm))
	ldps = []
	for wcc in wccs
		if length(wcc) >= min_size
			wcc_ldps, wcc_allele = optimize_paths(chrm[wcc,:], g2[wcc], W2[wcc,wcc], min_size, max_strands)
			allele[wcc] .= wcc_allele
			for wcc_ldp in wcc_ldps
				push!(ldps, wcc_ldp)
			end
		end
	end
	ldps, allele
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

    max_locus = maximum(chr[!,"g"])
    for (locus, nbrs) in enumerate(nbrs)
        # add outgoing edges to genomically down stream spatial neighbors
        for nbr in nbrs
            if chr[nbr,"g"] > chr[locus,"g"]
                #genomic_weight = max_locus - (chr[nbr,"g"] - chr[locus,"g"])
				genomic_weight = 1/abs(chr[nbr,"g"] - chr[locus,"g"])
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
    #A2 = spzeros(nnodes+2, nnodes+2)
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

#function optimize_paths(chrm, g :: DiGraph, g2:: DiGraph, W :: SparseMatrixCSC, min_size :: Int64, max_strands :: Int64)
function optimize_paths(chrm, g:: DiGraph, W :: SparseMatrixCSC, min_size :: Int64, max_strands :: Int64)
	nbranches = 2
	n_locus_nodes = nv(g)
	n_nodes = n_locus_nodes + 2*nbranches*max_strands
	src_nodes = Array((n_locus_nodes+1):(n_locus_nodes + nbranches*max_strands))
	dst_nodes = Array((src_nodes[end]+1):(src_nodes[end] + nbranches*max_strands))
	add_vertices!(g, 2*nbranches*max_strands)
	imag_edges = []
	g_locus_only = copy(g)
	for i in 1:n_locus_nodes
		for src in src_nodes
			@assert add_edge!(g, src, i)
			push!(imag_edges, (src, i))
		end
		for dst in dst_nodes
			@assert add_edge!(g, i, dst)
			push!(imag_edges, (i, dst))
		end
	end

	model = Model(GLPK.Optimizer)
	#A2 = Graphs.LinAlg.adjacency_matrix(g)
	#rows, cols, vals = findnz(A2)

	g_edges = Tuple.(collect(Graphs.edges(g)))
	g_locus_edges = filter(e -> e[1] <= n_locus_nodes && e[2] <= n_locus_nodes, g_edges)

	@variable(model, x[g_edges], Bin, container=SparseAxisArray)
	@variable(model, allele[1:(n_locus_nodes+2*nbranches*max_strands), 1:max_strands], Bin)
	#@variable(model, allele[1:(n_locus_nodes + 4*max_strands + 1), 1:max_strands], Bin)
	@variable(model, edge_allele[1:length(g_locus_edges),1:max_strands], Bin)

	# We set constraints to ensure that 

	#each node can have no more than one allele
	@constraint(model, sum(allele, dims=2) .<= 1)

	#set the src and dst node alleles
	for _allele in 1:max_strands
		for i in 1:nbranches, nodes in [src_nodes, dst_nodes]
		#for i in 1:1, nodes in [src_nodes, dst_nodes]
			#@constraint(model, allele[nodes[2*(_allele-1)+i], _allele] <= 1)
			@constraint(model, 1 <= allele[nodes[nbranches*(_allele-1)+i], _allele])
		end
	end

	#src and dst can have at most 1 edge
	#@constraint(model, [src in src_nodes], sum(x[(src,nbr)] for nbr in outneighbors(g, src)) <= 1)
	for src in src_nodes
		@constraint(model, sum(x[(src,nbr)] for nbr in outneighbors(g, src)) <= 1)
	end
	@constraint(model, [dst in dst_nodes], sum(x[(nbr,dst)] for nbr in inneighbors(g, dst)) <= 1)

	# each node must have the same allele as its parents
	#@constraint(model, [i = 1:n_nodes, nbr in inneighbors(g, i)], (x[(nbr,i)]-1) .<= allele[i,:] - allele[nbr,:])
	#@constraint(model, [i = 1:n_nodes, nbr in inneighbors(g, i)], allele[i,:] - allele[nbr,:] .<= (1-x[(nbr,i)]))

	# each node must have the same allele as its children
	for i in 1:n_locus_nodes, nbr in outneighbors(g, i)
	#for i in 1:n_nodes, nbr in outneighbors(g, i)
		@constraint(model, (x[(i,nbr)]-1) .<= allele[i,:] - allele[nbr,:])
		@constraint(model, allele[i,:] - allele[nbr,:] .<= (1-x[(i,nbr)]))
	end
	for i in src_nodes, nbr in outneighbors(g, i)
		@constraint(model, (x[(i,nbr)]-1) .<= allele[i,:] - allele[nbr,:])
		@constraint(model, allele[i,:] - allele[nbr,:] .<= (1-x[(i,nbr)]))
	end

	# if a locus node has no parents, it has no allele
	@constraint(model, [i = 1:n_locus_nodes], sum(allele[i,:]) <= sum(x[(nbr, i)] for nbr in inneighbors(g, i)))
	
	# an edge has the same allele as its parent node if the edge is active, no allele otherwise
 	for (i, edge) in enumerate(g_locus_edges)
 		@constraint(model, allele[edge[1],:] .+ (x[edge] - 1) .<= edge_allele[i,:])
 		@constraint(model, edge_allele[i,:] .<= allele[edge[1],:] .+ (1 - x[edge]))
		@constraint(model, edge_allele[i, :] .<= x[edge])
 	end

	# locus nodes have the same number of incoming and outgoing edges
	@constraint(model, [i = 1:n_locus_nodes], sum(x[(nbr,i)] for nbr in inneighbors(g, i)) == sum(x[(i,nbr)] for nbr in outneighbors(g, i)))

	# locus nodes may have no more than 2 incoming or outgoing edges
	@constraint(model, [i = 1:n_locus_nodes], sum(x[(nbr,i)] for nbr in inneighbors(g, i)) <= nbranches)
	#@constraint(model, [i = 1:n_locus_nodes], sum(x[(i,nbr)] for nbr in outneighbors(g, i)) <= 2)

	# no locus node can be connected to more than three other locus nodes
	@constraint(model, [i = 1:n_locus_nodes], sum(x[(nbr,i)] for nbr in inneighbors(g_locus_only, i)) + sum(x[(i,nbr)] for nbr in outneighbors(g_locus_only, i)) <= 3)

	# each src node can have at most one outgoing edge
	for src in src_nodes
		@constraint(model, sum(x[(src, nbr)] for nbr in 1:n_locus_nodes) <= 1)
	end

	# each dst node can have at most one incoming edge 
	@constraint(model, [dst in dst_nodes], sum(x[(nbr, dst)] for nbr in 1:n_locus_nodes) <= 1)

	# the two src nodes for an allele cannot start separate strands
	@constraint(model, sum(allele[1:n_locus_nodes,:], dims=1) .<= sum(edge_allele, dims=1) .+ 1)

	@objective(model, Max, sum(x[e]*W[e...] for e in g_locus_edges) - sum(x[e] for e in imag_edges))

	optimize!(model)

	evals = value.(x)
	
	gres = DiGraph(n_locus_nodes)
	for e in g_edges
		#evals[e] == 1 ? add_edge!(gres, e[1], e[2]) : nothing
		evals[e] == 1 ? add_edge!(gres, e) : nothing
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

function cluster_chromosomes_DBSCAN(chrms, radius, min_neighbors, min_size)
	points = Array([chrms.x chrms.y chrms.z]')
	dbr = dbscan(points, radius, min_neighbors=min_neighbors, min_cluster_size=min_size)
	dbscan_clusters = [sort(vcat(dbc.core_indices,  dbc.boundary_indices)) for dbc in dbr]
	return dbscan_clusters
end

#function reconcile_LDP_DBSCAN(ldp_allele, dbscan_clusters)
function compare_LDP_DBSCAN(ldps, dbscan_clusters, chrm, max_strands, ldp_overlap_thresh)

	#LDPs = [Array(1:nrow(chrms))[chrms.allele .== p] for p in 1:maximum(ldp_allele)]

	# find intersections of clusters by dbscan and ldp
	intersections = [ldp ∩ dbc for ldp in ldps, dbc in dbscan_clusters]

	#find the size of the intersections of clusters found by dbscan and ldp
	intersection_sizes = length.(intersections)

	accepted = zeros(Bool, size(intersections)...)

	#check whether the cluster agreement passes the threshold
	for ldp_ind in 1:max_strands
	    accepted[ldp_ind, :] = (intersection_sizes[ldp_ind,:]./length(ldps[ldp_ind]) .> ldp_overlap_thresh)
	end

	# find the number of DBSCAN clusters that each LDP cluster overlaps with
	n_ldp_dbc_overlaps = sum(intersection_sizes .> 0, dims=2)
	# find the number of LDP clusters that each DBSCAN cluster overlaps with
	n_dbc_lpc_overlaps = sum(intersection_sizes .> 0, dims=1)

	# Every LDP cluster must match two a DBSCAN cluster, and every DBSCAN cluster must match to 1 or 0 LDP cluster
	# for us to use the consensus of the two methods. Otherwise, we use LDP only.
	if (all(n_ldp_dbc_overlaps .== 1) && all(n_dbc_lpc_overlaps .<= 1))
	    #use union of ldp and its overlapping dbc
		ldp, dbc, val = findnz(sparse(accepted))
	    clusters = [sort(ldps[ldp[i]] ∪ dbscan_clusters[dbc[i]]) for i in 1:length(ldp)]
		allele = get_allele_col(chrm, clusters)
	else
		#use ldps
		allele = get_allele_col(chrm, ldps)
	end
	return allele
end
