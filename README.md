# DNASeqFISHChromosomeAssignment

This is a julia package for assigning loci in DNASeqFISH to chromosomes. The package uses two algorithms to make the assignments: [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) and Longest Disjoint paths which we introduce in our paper.

The package provides the function <code> assign_chromsomes </code>. Its use is documented in its docstring in [assign_functions.jl](https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment/blob/main/src/assignment_functions.jl).
