# DNASeqFISHChromosomeAssignment.jl

This is a julia package for assigning loci in DNASeqFISH to chromosomes. The package uses two algorithms to make the assignments: [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) and Longest Disjoint Paths which we introduce in "Separation of Homologous Chromosomes" subsection of the methods section in our preprint, [High-resolution spatial multi-omics reveals cell-type specific nuclear compartments](https://www.biorxiv.org/content/10.1101/2023.05.07.539762v1.abstract). You can find the main github repository for the data processing and analysis for the preprint [here](https://github.com/CaiGroup/dna-seqfish-plus-multi-omics/tree/main).

The package provides the function <code>assign_chromsomes</code>. Its use is documented in its docstring in [assign_functions.jl](https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment/blob/main/src/assignment_functions.jl).

For more information on how to use the functions provided in this package, see the [documentation](https://caigroup.github.io/DNASeqFISHChromosomeAssignment.jl/). The [test folder](https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment.jl/tree/main/test) contains a [script](https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment.jl/blob/main/test/test_e14_data.jl) that shows an example of how to use the package with sample example data.

If you use this package in a journal article, please cite our [preprint](https://www.biorxiv.org/content/10.1101/2023.05.07.539762v1.abstract).


# Installation

from the [julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) type
```
julia> ] add https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment.jl/
```
Or if in the root directory of the repostiory
```
julia> ] add .
```
