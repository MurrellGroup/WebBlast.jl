#Crude parsing of the BLAST XML. I could make types for Interations/Queries, Hits, and HSPs, but...

"""
    blast_xml_to_dict_arr(xmldoc::EzXML.Document)

```

[Note: if you CBA to figure out how to traverse all this, just call ```flatten_to_dataframe(blast_arr)```]

Takes an XML document, returned from BLAST and pushed through EzXML, and converts this to an Array of dictionaries.

The base array has one element per query (with one query, this will be an array of length 1).

Each Query dictionary will have a key "Query_hit_array", which will return an array of hit dictionaries.

Each hit can have more than one "High Similarity Pair" (or "HSP"), so the hit dictionary has a key "Hit_hsps" which itself is an array of HSPs.

```
"""
function blast_xml_to_dict_arr(xmldoc::EzXML.Document)
    blast_arr = []
    iteration_nodes = findall("//Iteration", xmldoc.root)
    for iter_n in iteration_nodes
        iter_dict = Dict()
        for iter_els in elements(iter_n)
            if iter_els.name == "Iteration_hits"
                iter_dict["Query_hit_array"] = []
                for h in elements(iter_els)
                    hit_dict = Dict()
                    for e in elements(h)
                        if e.name == "Hit_hsps"
                            hit_dict["Hit_hsps"] = []
                            for hsp in elements(e)
                                d = Dict([el.name => el.content for el in elements(hsp)])
                                push!(hit_dict["Hit_hsps"],d)
                            end
                        else
                            hit_dict[e.name] = e.content
                        end
                    end
                    push!(iter_dict["Query_hit_array"],hit_dict)
                end
            else        
                if !(iter_els.name in ["Iteration_query-ID","Iteration_stat"])
                    iter_dict[replace(replace(iter_els.name,"Iteration" => "Query"),"_query" => "")] = iter_els.content
                end
            end
        end
        push!(blast_arr,iter_dict)
    end
    return blast_arr
end

export flatten_to_dataframe

"""
    flatten_to_dataframe(blast_arr)

Takes an output array of dictionaries (ie. something returned by a ```WebBLAST(query)``` call) and flattens into a DataFrame.
"""
function flatten_to_dataframe(blast_arr)
    #The fields to keep are the things in here, minus the array fields
    #println(union(vcat(collect.(keys.(blast_arr))...)))
    #println(union(vcat(collect.(keys.(blast_arr[1]["Query_hit_array"]))...)))
    #println(union(vcat(collect.(keys.(blast_arr[1]["Query_hit_array"][1]["Hit_hsps"]))...)))
    query_fields = ["Query-def", "Query_iter-num", "Query-len"]
    hit_fields = ["Hit_accession", "Hit_def", "Hit_len", "Hit_num", "Hit_id"]
    HSP_fields = ["Hsp_query-from", "Hsp_align-len", "Hsp_hseq", "Hsp_hit-from", "Hsp_score", "Hsp_positive", "Hsp_qseq", "Hsp_query-to", "Hsp_query-frame", "Hsp_midline", "Hsp_evalue", "Hsp_hit-to", "Hsp_identity", "Hsp_hit-frame", "Hsp_gaps", "Hsp_bit-score", "Hsp_num"];
    flat_fields = vcat(query_fields,hit_fields,HSP_fields);

    field_arr_of_arr = [[] for i in 1:length(flat_fields)]
    for query in blast_arr
        for hit in query["Query_hit_array"]
            for HSP in hit["Hit_hsps"]
                ind = 1
                for f in query_fields
                    push!(field_arr_of_arr[ind],get(query,f,""))
                    ind += 1
                end
                for f in hit_fields
                    push!(field_arr_of_arr[ind],get(hit,f,""))
                    ind += 1
                end
                for f in HSP_fields
                    push!(field_arr_of_arr[ind],get(HSP,f,""))
                    ind += 1
                end
            end
        end
    end

    df = DataFrame()
    for i in 1:length(flat_fields)
        df[!,flat_fields[i]] = field_arr_of_arr[i]
    end
    return df
end

export WebBLAST

"""
    WebBLAST(query;
        num_hits = 100,
        database = "nt",
        program = "blastn",
        verbosity = 1, 
        option_string = "",
        save_XML_path = nothing,
        blast_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?")

```

NOTE: This relies on querying a public web service, using an API that they list as "deprecated". Please do not count on this working in the long run, and please don't abuse it as that might get it shut down faster.

WHAT YOU PUT IN:

query must be either a String, or an array of String (if searching multiple sequences). These are expected to be nucleotides, or amino acids, but we aren't bothering with specific types for these. Just make sure they match the database and program.

save_XML_path will export the returned XML - you probably don't need this.

num_hits controls the maximum number of hits that get returned.

database and program are common blast options you might want to modify. Eg. program = "blastp", database = "pdb".

option_string is any string you want to get injected into the BLAST Put call. Must be of the form "&MATRIX=PAM230&NUCL_REWARD=2" etc 

blast_URL is a string that you can use to point to a different blast service, if you have access to one.

WHAT GETS RETURNED:

[Note: if you CBA to figure out how to traverse all this, just call ```flatten_to_dataframe(WebBLAST(query))``` instead.]

Takes an XML document, returned from BLAST and pushed through EzXML, and converts this to an Array of dictionaries.

The base array has one element per query (with one query, this will be an array of length 1).

Each Query dictionary will have a key "Query_hit_array", which will return an array of hit dictionaries.

Each hit can have more than one "High Similarity Pair" (or "HSP"), so the hit dictionary has a key "Hit_hsps" which itself is an array of HSPs.
```

"""
function WebBLAST(query;
        num_hits = 100,
        database = "nt",
        program = "blastn",
        verbosity = 1, 
        option_string = "",
        save_XML_path = nothing,
        blast_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?")
    #This could be done with dispatch, but then I'd have to have all the options repeated.
    if typeof(query) == Vector{String}
        q = ""
        for (i,qu) in enumerate(query)
            #Note URL encoding of > (%3E) and \n (%0A)
            q *= "%3Equery_$(i)%0A"*qu*"%0A"
        end
    elseif typeof(query) == String
        q = query
    else
        @error "query must either be a String or a vector of Strings"
    end
    #Unclear that the "XML" FORMAT_TYPE is working. But no worries for now.
    qstr = "CMD=Put&QUERY="*q*"&PROGRAM=$(program)&DATABASE=$(database)&FORMAT_TYPE=XML&HITLIST_SIZE=$(num_hits)"*option_string;    
    resp = HTTP.request("POST", blast_URL,
             ["Content-Type" => "application/x-www-form-urlencoded"],
             qstr);
    
    #resp = HTTP.post(blast_URL*qstr);
    respstr = String(HTTP.body(resp));
    RID = ""
    RTOE = 0
    try
        RID = match(r"RID \= \w*\n",respstr).match[7:end-1]
        RTOE = max(20,parse(Int64,match(r"RTOE \= \w*\n",respstr).match[8:end-1]))
    catch
        @error "Problem during initial call. Something is wrong with your query/queries, or they are incompatible with the program/db you're trying to use."
        return nothing
    end
    verbosity > 0 && println("Waiting for $(RTOE) seconds")
    sleep(RTOE)
    http_str = blast_URL*"CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$(RID)"
    getresp = HTTP.get(http_str);
    getstr = String(HTTP.body(getresp));
    
    #loop that waits for BLAST to finish
    while occursin("Status=WAITING",getstr)
        verbosity > 0 && println("Waiting for result...")
        sleep(20)
        getresp = HTTP.get(http_str);
        getstr = String(HTTP.body(getresp))
    end
    
    if occursin("Status=READY",getstr) && occursin("ThereAreHits=yes",getstr)
        verbosity > 0 && println("Retrieving hits.")
        xmlblast = HTTP.get(blast_URL*"CMD=Get&FORMAT_TYPE=XML&RID=$(RID)")
        xmlstr = String(HTTP.body(xmlblast))
        if typeof(save_XML_path) == String
            open(save_XML_path, "w") do io
               write(io, xmlstr)
            end
        end
        xmldoc = EzXML.parsexml(xmlstr)
        return blast_xml_to_dict_arr(xmldoc)
    else
        if verbosity > 1
            @warn "ERROR - BLAST failed. Returning latest HTTP.get request string instead of BLAST results, just in case you can figure out what went wrong."
            return getstr
        else
            @error "BLAST failed. If verbosity is set higher than 1, the latest HTTP.get request string will be returned to allow debugging."
        end
    end
end