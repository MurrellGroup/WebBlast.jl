using WebBlast, EzXML, DataFrames
using Test

@testset "WebBlast.jl" begin

    read_xml = String(read("test.xml"))
    xmldoc = EzXML.parsexml(read_xml)
    blast_arr = WebBlast.blast_xml_to_dict_arr(xmldoc)
    @test length(blast_arr) == 2
    df = flatten_to_dataframe(blast_arr)
    @test size(df) == (4, 25)
end
