# WebBlast.jl

This package helps you call BLAST (https://blast.ncbi.nlm.nih.gov/) from Julia. It is extremely simple, but tries to remain flexible, as long as you're happy to put in a little extra work.

The main functions you'll use is:

```julia
    WebBLAST(query;
        num_hits = 100,
        database = "nt",
        program = "blastn")
```

NOTE: This relies on querying a public web service, using an API that they list as "deprecated". Please do not count on this working in the long run, and please don't abuse it as that might get it shut down faster.

```query``` must be either a String, or an array of String (if searching multiple sequences). These are expected to be nucleotide or amino acid sequences. Just make sure they (nuc or AA) match the database and program! For example, if you pass in amino acid sequences, you can't use the "nt" database, but database = "pdb" should work.

## Installation:
Until the package is registered, please use:
```julia
using Pkg
Pkg.add(PackageSpec(name="WebBlast", url = "https://github.com/MurrellGroup/WebBlast.jl.git"))
```

## Example use:
```julia
using WebBlast

queries = ["MKTLLLTIGFSLIAILQAQDTPALGKDTVAVSGKWYLKAMTADQEVPE","MRALLLAIGLGLVAALQAQEFPAVGQPLQDLLGRWYLKAMTSDPEIPG"];
blast_array = WebBLAST(queries, program = "blastp", database = "pdb", num_hits = 3);
```
with output (note the delays, while the BLAST server computes):
```
Waiting for 35 seconds
Waiting for result...
Waiting for result...
Retrieving hits.
```
and returns ```blast_array```, which is an array of dictionaries. That array has one element per query (with one query, this will be an array of length 1). Each Query dictionary will have a key "Query_hit_array", which will return an array of hit dictionaries. Each hit can have more than one "High Similarity Pair" (or "HSP"), so the hit dictionary has a key "Hit_hsps" which itself is an array of HSPs. This can be traversed to get whatever you like. For example, here is a code snippet that runs through all HSPs for all hits for all queries and prints out some key properties of the hits:
```julia
for query in blast_array
    for hit in query["Query_hit_array"]
        for HSP in hit["Hit_hsps"]
            println("Query: ", query["Query-def"])
            println("Accession: ", hit["Hit_accession"])
            println("Identity: ", HSP["Hsp_identity"])
            println("Query: ", HSP["Hsp_qseq"])
            println("Hit:   ", HSP["Hsp_hseq"])
        end
    end
end
```
Which gives the following output:
```
Query: query_1
Accession: 5T43_A
Identity: 14
Query: VSGKWYLKAMTADQEVPE
Hit:   VSGTWYLKAMTVDREFPE
Query: query_1
Accession: 1XKI_A
Identity: 14
Query: VSGKWYLKAMTADQEVPE
Hit:   VSGTWYLKAMTVDREFPE
Query: query_1
Accession: 4RUN_A
Identity: 9
Query: VSGKWYLKAMTADQEVPE
Hit:   ITGTWYVKAMVVDKDFPE
Query: query_2
Accession: 5T43_A
Identity: 14
Query: AVGQPLQDLLGRWYLKAMTSDPEIP
Hit:   ASDEEIQDVSGTWYLKAMTVDREFP
Query: query_2
Accession: 1XKI_A
Identity: 14
Query: AVGQPLQDLLGRWYLKAMTSDPEIP
Hit:   ASDEEIQDVSGTWYLKAMTVDREFP
Query: query_2
Accession: 4QAF_A
Identity: 12
Query: AVGQPLQDLLGRWYLKAMTSD
Hit:   ASDEEIQDVSGTWYLKAMTVD
```

However, if you don't want to both with all that, you can instead use:
```julia
df = flatten_to_dataframe(blast_array)
```
which will give you a DataFrame with all of the information per HSP, hit, query, etc.
