# DNASeqFISHChromosomeAssignment.jl

This is a julia package for assigning loci in DNASeqFISH to chromosomes. The package uses two algorithms to make the assignments: [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) and Longest Disjoint paths which we introduce in our paper, [High-resolution spatial multi-omics reveals cell-type specific nuclear compartments](https://www.biorxiv.org/content/10.1101/2023.05.07.539762v1.abstract). You can find the main page for the data processing and analysis for the paper [here](https://github.com/CaiGroup/dna-seqfish-plus-multi-omics/tree/main).


# Installation

from the julia REPL type
```
julia> ] add https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment.jl/
```
Or if in the root directory of the repostiory
```
julia> ] add .
```
The package provides the function <code>assign_chromsomes</code>. Its use is documented in its docstring in [assign_functions.jl](https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment/blob/main/src/assignment_functions.jl).

For more information, see the [documentation](https://caigroup.github.io/DNASeqFISHChromosomeAssignment.jl/).

If you use this package in a journal article, please cite our [preprint](https://www.biorxiv.org/content/10.1101/2023.05.07.539762v1.abstract).
