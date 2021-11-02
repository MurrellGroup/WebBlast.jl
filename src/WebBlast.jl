#Blegh. Despite reasonable prior googling, I somehow missed this:
#https://github.com/hng/BiomolecularStructures.jl/blob/master/src/WebBLAST/ncbi_blast.jl
#Some bits of that are nice, and some are inconvenient.
#Their restrictive assumptions are:
#the PDB is always the target.
#Sequences are always amino acids.
#You only ever search with one query
#Still, it would have been much better to start with that and adapt...
#Also, now the name WebBLAST is taken. Should I still call the package WebBlast.jl? Ja, why not.

module WebBlast
    using HTTP, EzXML, DataFrames
    include("functions.jl")
end
