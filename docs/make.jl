#push!(LOAD_PATH,"../src/")
using Documenter
using DNASeqFISHChromosomeAssignment

makedocs(
    sitename = "DNASeqFISHChromosomeAssignment",
    format = Documenter.HTML(),
    modules = [DNASeqFISHChromosomeAssignment]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
