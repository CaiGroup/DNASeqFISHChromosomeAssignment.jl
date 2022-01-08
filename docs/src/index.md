# DNASeqFISHChromosomeAssignment.jl

## Introduction

This package implements the Longest Disjoint Paths algorithm for identifying chromosome strands in DNA SeqFISH+ data.
The intuition behind this method DNA SeqFISH+ is like imaging beads on an invisible string, where the loci we measure are the beads and the genomic DNA is the string.
We know the spatial and genomic coordinates of each locus we find, so we can infer that DNA strands connect loci that are close both spatially and genomically.

However, we do not try to model the branching nature of replicating chromosomes. To correctly classify replicated loci on chromosome pairs that are spatially distinct, we combine
the LDP chromosome identification method with DBSCAN.

## Installation



## API Reference

```@docs
assign_chromosomes
```