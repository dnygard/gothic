# full run of the pipeline, from sam files to final results dataframe plus figures
# files that need to be in filedir:
# - map(HUMAN_9606_idmapping_selected.tsv)
# - gencode(gencode_pcgenes.csv)
# - go obo(go-basic.obo)
# - GothicAPSP.so
# - Ensembl2Reactome_All_Levels.txt

from gothic import *
import sys
from gprofiler import GProfiler


# performs spectral clustering on graph
def do_spectral_clustering(graph, outfile, method="OPTICS"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    # get graph spectrum for clustering
    adjacency_matrix = adjacency(g, weight=g.ep.weight).toarray()  # Compute the weighted adjacency matrix
    degree_matrix = np.diag(adjacency_matrix.sum(axis=0))  # Compute the degree matrix for weighted graph
    laplacian_matrix = degree_matrix - adjacency_matrix  # Compute the Laplacian matrix
    eigenvalues, eigenvectors = eigsh(laplacian_matrix, k=3, which='SM')  # Eigenvalue Decomposition

    # do OPTICS Clustering
    clust = OPTICS(min_samples=3, xi=0.05).fit(eigenvectors)  # do clustering
    # cluster_labels = clust.labels_[clust.ordering_]
    cluster_labels = clust.labels_

    clusternodedict = {}
    node = 0
    for i in [int(x) for x in cluster_labels]:
        if i == -1:
            node = node + 1
            continue
        else:
            if i in clusternodedict.keys():
                clusternodedict[i].append(node)
            else:
                clusternodedict[i] = [node]
            node = node + 1

    with open(outfile, "w") as f:
        for k in clusternodedict.keys():
            f.write(",".join([str(x) for x in clusternodedict[k]]) + "\n")


# performs Markov clustering on graph
def do_markov_clustering(graph, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    matrix = adjacency(g, weight=g.ep.raw_weight).toarray()

    result = mcl.run_mcl(matrix)  # run MCL with default parameters
    clusters = mcl.get_clusters(result)  # get clusters

    for line in clusters:
        clusterhist.append(len(line))
        writeline = str(line).split("(")[-1].split(")")[0] + "\n"
        f.write(str(writeline))


# gets enrichment of GO terms for all clusters in cluster file from
def get_go_enrichments(graph, clustersfile, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    # read in cluster file
    with open(clustersfile, "r") as f:
        resultsDFlist = []
        for line in f:  # each line is a cluster
            splitlist = [int(x) for x in line.split(",")[:-1]]  # get list of nodes in cluster
            genelist = []
            for node in splitlist:
                genelist.extend(g.vp.genes[node])

            # for cluster in file
            gp = GProfiler(return_dataframe=True)
            if genelist:
                newDF = gp.profile(organism='hsapiens', query=genelist)
                newDF.to_csv(outfile, mode='a')


def gothic_full_clustering_run(runname, sam1=None, sam2=None, filedir=".", binsize=80000, vmax=10000, ncores=1,
                    step=0, endstep=12, binfile=None, saveadjlist=False, filterlevel=0.5):
    start = timeit.default_timer()

    if binfile:
        bindict = bins_to_dict(binfile)
        print("Bin dict created from file...")
    else:
        bindict = None

    # create graph
    if step <= 1 <= endstep:
        print("Creating Graph...")
        adjlist = sam_adjacencies(sam1, sam2, m=binsize, bins=bindict, verbose=True)
        print("Adjacency list created...")
        if saveadjlist:  # write adjacency list to file if saveadjlist option is true
            print("Writing adjacency list to file...")
            adjlistoutfile = runname + "_ADJLIST.tsv"
            save_adjlist(adjlist, adjlistoutfile)

    # convert adjlist to graphml
    if step <= 2 <= endstep:
        print("Converting adjlist to .graphml...")
        graphname = runname + "_noannotation.graphml"
        if saveadjlist:
            adjlist = []  # lil hack to overwrite the adjlist or else create one if there is none (ie start w this step)
            del adjlist  # clear up memory from filled adjlist dictionary object
            adjlistoutfile = runname + "_ADJLIST.tsv"
            adjlist_file_to_graphml(adjlistoutfile, outfile=graphname)
        else:
            hic_adjlist_to_graphml(adjlist, fileprefix=graphname)
            del adjlist  # clear up memory from filled adjlist dictionary object
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(GraphStep): ', stop - start)  # Timer Block

    # annotate graph
    if step <= 3 <= endstep:
        print("Annotating graph...")
        graphname = runname + "_noannotation.graphml"
        mymap = filedir + "/HUMAN_9606_idmapping_selected.tsv"
        gencode = filedir + "/gencode_pcgenes.csv"
        goobo = filedir + "/go-basic.obo"
        reactomemap = filedir + "/Ensembl2Reactome_All_Levels.txt"
        newgraphname = runname + "_incompleteannotated.graphml"
        genes_go_annotate(graphname, mapfile=mymap, gencodefile=gencode, m=binsize,
                          binfile=binfile, outfile=newgraphname, go_obo=goobo)
        oldgraphname = newgraphname
        newgraphname = runname + "_annotated.graphml"
        reactome_annotate(oldgraphname, mapfile=reactomemap, outfile=newgraphname)

        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(AnnotateStep): ', stop - start)  # Timer Block

    # normalize, adjust weights, filter, extract largest component
    if step <= 4 <= endstep:
        print("Normalizing...")
        graphname = runname + "_incompleteannotated.graphml"
        newgraphname = runname + "_normalized.graphml"
        cooname = runname + "prenormCOO.txt"
        normcooname = runname + "normedCOO.txt"
        prenorm_filter(graphname, graphname)  #filter edges with only one contact before normalizing
        make_graph_coo(graphname, outfile=cooname)
        ice_balance(cooname, normcooname)
        update_weights_from_coo(graphname, normcooname, newgraphname)

        print("filtering edges...")
        graphname = newgraphname
        fpstring = "Top" + str(int(filterlevel * 100)) + "p"
        newgraphname = runname + "_" + fpstring + "filtered.gt"
        create_top_edges_graph(graphname, threshold=filterlevel, outfile=newgraphname)
        print('Global Runtime(NormalizeStep): ', stop - start)  # Timer Block

    # plot edge weights and make graph report
    if step <= 5 <= endstep:
        graphname = newgraphname
        weightplotname = runname + "_edgeWeightPlot.png"
        generate_graph_report(graphname, outfile=runname + "graphReport.txt")
        plot_edge_weight_distribution(graphname, outfile=weightplotname)

    # do spectral optics and mcl clustering
    if step <= 6 <= endstep:
        graphname = newgraphname
        fpstring = "Top" + str(int(filterlevel * 100)) + "p"
        spectralresultsfile = runname + "_" + fpstring + "SpectralClusters.csv"
        mclresultsfile = runname + "_MCLclusters.csv"

        # do clustering
        do_markov_clustering(graphname, mclresultsfile)
        do_spectral_clustering(graphname, spectralresultsfile)

        # do clustering enrichment
        mcl_plotfile = runname + "_" + fpstring + "MCLclusters.csv"
        spectral_plotfile = runname + "_" + fpstring + "SpectralClusters.csv"

        get_go_enrichments(graphname, mclresultsfile, mcl_plotfile)
        get_go_enrichments(graphname, spectralresultsfile, spectral_plotfile)


if __name__ == "__main__":
    runname = sys.argv[1] + "250kbBins"
    samone = sys.argv[1] + "_1.sam"
    samtwo = sys.argv[1] + "_2.sam"

    gothic_full_clustering_run(runname, sam1=samone, sam2=samtwo, filterlevel=0.2, saveadjlist=True, binsize=250000)

