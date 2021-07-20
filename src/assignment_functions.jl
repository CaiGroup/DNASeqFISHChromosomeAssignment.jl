using DataFrames
using JuMP
using NearestNeighbors
using LightGraphs
using SparseArrays
using GLPK
using Clustering

"""
assign_chromosomes(pnts :: DataFrame,
	 			   search_radius :: Real,
				   sigma :: Real,
				   max_strands :: Int64,
				   min_size :: Int64 = 2)

Assigns genomic loci from DNA seqFISH data to indiviual chromosomes. Builds
a directed acyclic graph where every locus has an out going edge to genomically
down stream loci within the search radius. Edges are weighted by:

edge_weight = (max_genomic_coord - (parent_genomic coord - child_genomic_coord))*
			exp(-(||x_parent - x_child||^2)/(2σ^2))

Searches for up to max_strands disjoint paths that maximize the sum of
the weights of all edges in the paths

Takes a DataFrame (pnts) of decoded DNA seqFISH data where each column describes a locus.
Must include columns:
	"fov" : the field of view in which locus was imaged
	"cellID" : the cell in which the locus was imaged
	"chrom" : chromosome on which locus is located
	"x" : x spatial coordinate of locus
	"y" : y spatial coordinate of locus
	"z" : z spatial coordinate of locus
	"pos" : genomic coordinate of locus

The search radius gives how far in space to search for neighboring loci on the same
chromosome strand for each locus.

The sigma value defines the weight function for how much to penalize grouping adjacent
genomic loci onto the the same chromosome based on their spatial distance from each other

The minimum size is the minimum acceptable number of loci on a strand.

returns a copy of the input points sorted by fov, cell_ID, chrom, and pos with a
new column: allele.

allele = -1 means that a locus was not assigned to a chromosome.
allele > 0 indicates the chromsome number that the locus was assigned to.
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
	sort!(pnts, [:fov, :cellID, :chrom, :pos])
	chrms = groupby(pnts,[:fov, :cellID, :chrom])
	assn_chrms(chrm) = assign_loci(chrm, search_radius, sigma, max_strands, min_size, dbscan, dbscan_min_pnts, ldp_overlap_thresh)
	res = transform(assn_chrms, chrms)
	#renamed_col = combine(res, Not(:x1), :x1 => :allele)
	return res #renamed_col
end

function assign_loci(chrm, r :: Real, sig :: Real, max_strands :: Int64, min_size :: Int64, dbscan, dbscan_min_pnts, ldp_overlap_thresh)
	#chrm = sort(_chrm, :pos)
	println("fov: ", chrm[1, "fov"], ", cell: ", chrm[1, "cellID"], ", ", chrm[1,"chrom"])
	ldps, ldp_allele = find_longest_disjoint_paths(chrm, r, sig, max_strands, min_size)
	return_df = DataFrame(Dict("ldp_allele"=>ldp_allele))
	if dbscan
		dbscan_clusters = cluster_chromosomes_DBSCAN(chrm, r, dbscan_min_pnts, min_size)
		return_df[!, "dpscan_allele"] = get_allele_col(chrm, dbscan_clusters)
		return_df[!, "final_allele"] = compare_LDP_DBSCAN(ldps, dbscan_clusters, chrm, max_strands,ldp_overlap_thresh)
	end
	return return_df
end

function find_longest_disjoint_paths(chrm, r :: Real, sig :: Real, max_strands :: Int64, min_size :: Int64)
	A, W, g = get_neighbors(chrm, r, sig)
	R, WR = get_trans_red_w(chrm, r, sig, A, W, g)
	g2, W2 = rule_out_edges(A, W, WR)
	return optimize_paths(chrm, g, g2, W2, min_size, max_strands)
end

function get_neighbors(chr, r, sig)
    nloci = size(chr)[1]
    chr = DataFrame(chr)
    sort!(chr, "pos")

    # get spatial neighbors
    coords = Array([chr.x chr.y chr.z]')
    tree = KDTree(coords)
    nbrs = inrange(tree, coords, r)

    # build neighbor adjacency graph weighted by genomic distance
    A = spzeros(Bool, nloci, nloci)
    W = spzeros(Float64, nloci, nloci)
    g = SimpleDiGraph(nloci)

    max_locus = maximum(chr[!,"pos"])
    for (locus, nbrs) in enumerate(nbrs)
        # add outgoing edges to genomically down stream spatial neighbors
        for nbr in nbrs
            if chr[nbr,"pos"] > chr[locus,"pos"]
                genomic_weight = max_locus - (chr[nbr,"pos"] - chr[locus,"pos"])
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

	model = Model(GLPK.Optimizer)
	#A2 = LightGraphs.LinAlg.adjacency_matrix(g2)
	#rows, cols, vals = findnz(A2)

	g2_edges = Tuple.(collect(LightGraphs.edges(g2)))
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

	intersections = [ldp ∩ dbc for ldp in ldps, dbc in dbscan_clusters]
	intersection_sizes = length.(intersections)

	accepted = zeros(Bool, size(intersections)...)

	for ldp_ind in 1:max_strands
	    accepted[ldp_ind, :] = (intersection_sizes[ldp_ind,:]./length(ldps[ldp_ind]) .> ldp_overlap_thresh)
	end

	nz_intersections = (intersection_sizes .> 0)
	n_ldp_dbc_overlaps = sum(intersection_sizes .> 0, dims=2)
	n_dbc_lpc_overlaps = sum(intersection_sizes .> 0, dims=1)

	if ~ (all(n_ldp_dbc_overlaps .== 1) && all(n_dbc_lpc_overlaps .<= 1))
	    #use ldps
		allele = get_allele_col(chrm, ldps)
	else
	    #use union of ldp and its overlapping dbc
		ldp, dbc, val = findnz(sparse(accepted))
	    clusters = [sort(ldps[ldp[i]] ∪ dbscan_clusters[dbc[i]]) for i in 1:length(ldp)]
		allele = get_allele_col(chrm, clusters)
	end
	return allele
end
