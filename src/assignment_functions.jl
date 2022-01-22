using DataFrames
using JuMP
using NearestNeighbors
using Graphs
using SparseArrays
using GLPK
using Clustering
using CPLEX

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
	- x : x spatial coordinate of locus
	- y : y spatial coordinate of locus
	- z : z spatial coordinate of locus
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
	 						r_dbscan :: Real,
							r_ldp :: Real, 
							sigma :: Real,
							min_size :: Int64 = 100,
							min_prop_unique :: Float64 = 0.9,
							dbscan_min_pnts :: Int64 = 10)

	@assert min_size > 1
	sort!(pnts, [:fov, :cellID, :chrom, :g])
	chrms = groupby(pnts,[:fov, :cellID, :chrom])
	assn_chrms(chrm) = assign_loci(chrm, r_dbscan, r_ldp, sigma, min_size, min_prop_unique, dbscan_min_pnts)
	res = transform(assn_chrms, chrms)
	#renamed_col = combine(res, Not(:x1), :x1 => :allele)
	return res #renamed_col
end

function assign_loci(chrm, r_dbscan :: Real, r_ldp :: Real, sig :: Real,min_size :: Int64, min_prop_unique, dbscan_min_pnts)
	#chrm = sort(_chrm, :g)
	println("fov: ", chrm[1, "fov"], ", cell: ", chrm[1, "cellID"], ", ", chrm[1,"chrom"])
	dbscan_clusters = cluster_chromosomes_DBSCAN(chrm, r_dbscan, dbscan_min_pnts, min_size)
	dbscan_allele = get_allele_col(chrm, dbscan_clusters)
	final_allele = copy(dbscan_allele)
	for (i, c) in enumerate(dbscan_clusters)
		println("prop unique: ", length(unique(chrm.g[c]))/length(chrm.g[c]))
		if length(unique(chrm.g[c]))/length(chrm.g[c]) < min_prop_unique
			"starting ldp"
			tangled_chrms = chrm[c, :]
			ldps, ldp_allele = find_longest_disjoint_paths(tangled_chrms, r_ldp, sig, 2, min_size)
			if length(ldps) == 2
				println("found 2 ldps")
				ldp_allele[ldps[1]] .= i
				ldp_allele[ldps[2]] .= (maximum(dbscan_allele) + 1)
				final_allele[c] .= ldp_allele
			end
		end
	end
	#ldps, ldp_allele = find_longest_disjoint_paths(chrm, r, sig, max_strands, min_size)
	#return_df = DataFrame(Dict("ldp_allele"=>ldp_allele))
	return_df = DataFrame(Dict("dbscan_allele"=>dbscan_allele, "final_allele"=>final_allele))
	return return_df
end

function find_longest_disjoint_paths(chrm, r :: Real, sig :: Real, max_strands :: Int64, min_size :: Int64)
	A, W, g = @time get_neighbors(chrm, r, sig)
	R, WR = @time get_trans_red_w(chrm, r, sig, A, W, g)
	g2, W2 = @time rule_out_edges(A, W, WR)
	return @time optimize_paths(chrm, g, g2, W2, min_size, max_strands)
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
    A2 = spzeros(nnodes+2, nnodes+2)
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

function optimize_paths(chrm, g :: DiGraph, g2:: DiGraph, W :: SparseMatrixCSC, min_size :: Int64, max_strands :: Int64)
	src = nv(g)+1
	dst = src+1
	for i in 1:nv(g2)
		add_edge!(g2, src, i)
		add_edge!(g2, i, dst)
	end

	#model = Model(GLPK.Optimizer)
	model = Model(CPLEX.Optimizer)
	#A2 = Graphs.LinAlg.adjacency_matrix(g2)
	#rows, cols, vals = findnz(A2)

	g2_edges = Tuple.(collect(Graphs.edges(g2)))
	g2_locus_edges = filter(e -> e[1] <= nv(g) && e[2] <= nv(g), g2_edges)

	@variable(model, x[g2_edges], Bin, container=SparseAxisArray)
	@constraint(model, sum(x[(src,nbr)] for nbr in outneighbors(g2, src)) <= max_strands)
	@constraint(model, sum(x[(nbr,dst)] for nbr in inneighbors(g2, dst)) <= max_strands)
	@constraint(model, [i = 1:nv(g)], sum(x[(nbr,i)] for nbr in inneighbors(g2, i)) == sum(x[(i,nbr)] for nbr in outneighbors(g2, i)))
	@constraint(model, [i = 1:nv(g)], sum(x[(nbr,i)] for nbr in inneighbors(g2, i)) <= 1)

	@objective(model, Max, sum(x[e]*W[e...] for e in g2_locus_edges))

	optimize!(model)

	evals = value.(x)
	gres = DiGraph(nv(g))
	for e in g2_edges
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
