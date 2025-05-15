from gothic import *

# a comprehensive script for automatically performing all gothic analyses from any step


def gothic_autorun(infile, runname, infile2=None, step=0, binsize=250000, binfile=None, filter="simple", filterlevel=1,
                   logging=True, verbose=False, filedir="."):

    # determine type of input file and go to appropriate step
    file_extension = infile.split(".")[-1]

    if file_extension == "fastq":
        print("FASTQ format not supported in autorun."
              + "Please align files in SAM format to Hg38.13 using an aligner such as Bowtie2 before using autorun.")

    elif file_extension == "sam":
        print("SAM file detected. Running graph creation from SAM...")
        # start from graph creation
        adjlist = sam_adjacencies(sam1, sam2, m=binsize, bins=bindict, verbose=True)
        print("Adjacency list created...")
        hic_adjlist_to_graphml(adjlist, fileprefix=runname + ".gt")
        del adjlist  # clear up memory from filled adjlist dictionary object
        # then proceed to normalization

    elif file_extension == "hic":
        print("hic file detected. Running graph creation from hic...")
        # start from graph creation
        tsvname = runname + "_HicEdgelist.tsv"
        edgelist_from_hic(infile, tsvname, binsize=binsize)
        graph_from_hic_edgelist(tsvfile, runname + ".gt")
        # then proceed to normalization

    elif file_extension == "pairs":
        print("pairs file detected. Running graph creation from pairs...")
        # start form graph creation
        create_graph_from_pairs(infile, runname + ".gt", binsize=binsize, verbose=verbose, outtsv=None)

    elif file_extension == "graphml":
        print("graphml file detected. Proceeding to specified step...")
        # check if there is a start step specified, otherwise proceed from graph normalization

    elif file_extension == "gt":
        print("gt file detected. Proceeding to specified step...")
        # check if there is a start step specified, otherwise proceed from graph normalization

    else:
        print("Unknown input file type. Please try again with a .sam, .hic, .pairs, .graphml, or .gt file")

    graphname = runname + ".gt"

    # create vector of strings graph property to track graph modifications/transformations
    g = load_graph(graphname)
    try:
        g.gp.mods
    except AttributeError:
        mods = g.new_gp("vector<string>")
        g.gp["mods"] = mods
    g.save(graphname)
    del g

    # filter
    if filter == "simple":
        print("filtering by number of contacts...")
        prenorm_filter(graphname, graphname, mincontacts=filterlevel+1)  # filters edges with contacts <= filterlevel
        # add modification to mods list graph property
        g = load_graph(graphname)
        try:
            g.gp.mods.append("Filtered: " + filter + " at filterlevel " + str(filterlevel))
        except AttributeError:
            mods = g.new_gp("vector<string>")
            g.gp["mods"] = mods
            g.gp.mods.append("Filtered: " + filter + " at filterlevel " + str(filterlevel))
        g.save(graphname)
        del g

    if filter == "top_percent":
        print("filtering to keep highest weighted edges by percentage...")
        create_top_edges_graph(graphname, threshold=filterlevel, outfile=graphname, weights_are_distances=False)
        # add modification to mods list graph property
        g = load_graph(graphname)
        try:
            g.gp.mods.append("Filtered: " + filter + " at filterlevel " + str(filterlevel))
        except AttributeError:
            mods = g.new_gp("vector<string>")
            g.gp["mods"] = mods
            g.gp.mods.append("Filtered: " + filter + " at filterlevel " + str(filterlevel))
        g.save(graphname)
        del g

    if filter == "bottom_percent":
        print("filtering to keep lowest weighted edges by percentage...")
        create_top_edges_graph(graphname, threshold=filterlevel, outfile=graphname,
                               weights_are_distances=True)
        # add modification to mods list graph property
        g = load_graph(graphname)
        try:
            g.gp.mods.append("Filtered: " + filter + " at filterlevel " + str(filterlevel))
        except AttributeError:
            mods = g.new_gp("vector<string>")
            g.gp["mods"] = mods
            g.gp.mods.append("Filtered: " + filter + " at filterlevel " + str(filterlevel))
        g.save(graphname)
        del g

    # normalize
    cooname = runname + "_COO.txt"
    normcooname = runname + "COO_Normalized.txt"
    make_graph_coo(graphname, outfile=cooname)
    ice_balance(cooname, normcooname)
    update_weights_from_coo(graphname, normcooname, graphname)
    # add modification to mods list graph property
    g = load_graph(graphname)
    try:
        g.gp.mods.append("Normalized with ICE")
    except AttributeError:
        mods = g.new_gp("vector<string>")
        g.gp["mods"] = mods
        g.gp.mods.append("Normalized with ICE")
    g.save(graphname)
    del g

    # annotate
    mymap = filedir + "/HUMAN_9606_idmapping_selected.tsv"
    gencode = filedir + "/gencode_pcgenes.csv"
    goobo = filedir + "/go-basic.obo"
    genes_go_annotate(graphname, mapfile=mymap, gencodefile=gencode, m=binsize,
                      binfile=binfile, outfile=graphname, go_obo=goobo)
    # add modification to mods list graph property
    g = load_graph(graphname)
    try:
        g.gp.mods.append("Annotated with genes")
        g.gp.mods.append("Annotated with goterms")
    except AttributeError:
        mods = g.new_gp("vector<string>")
        g.gp["mods"] = mods
        g.gp.mods.append("Annotated with genes")
        g.gp.mods.append("Annotated with goterms")
    g.save(graphname)
    del g

    # create graph report
    print("Generating graph report...")
    generate_graph_report(graphname, outfile=runname + "_graphReport.txt")

    # do TPD
    print("Running TPD analysis...")
    # calculate All Pairs Shortest Paths (APSP) and create distance matrix
    annotate_apsp(graphname, graphname)  # calculate all pairs shortest paths as graph annotations
    # add modification to mods list graph property
    g = load_graph(graphname)
    try:
        g.gp.mods.append("Annotated with shortest paths distances")
    except AttributeError:
        mods = g.new_gp("vector<string>")
        g.gp["mods"] = mods
        g.gp.mods.append("Annotated with shortest paths distances")
    g.save(graphname)
    del g
    matfilename = runname + "_distmatrix.tsv"
    distmatrix_from_annotations(graphname, outfile=matfilename)  # convert/format distances to distance matrix



    # do Spectral & MCL clustering

    #

if __name__ == "__main__":
    caselist = ["ENCFF929TYU", "ENCFF375UIA"]
    for case in caselist:
        gname = case + "250kb.gt"
        rname = case + "250kb"
        gothic_autorun(gname, rname)
