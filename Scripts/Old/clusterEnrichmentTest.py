from gothic import *
from gprofiler import GProfiler

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


if __name__ == "__main__":
    # for each case
    # for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
    for case in ["MEC1c"]:
        graphfile = case + "FixedBins500kb_edgeModifiedFilteredTop10p.gt"
        mcl_clustfile = case + "FixedBins500kbTop10pEdges_SpectralClusters.csv"
        mcl_outfile = case + "FixedBins500kbTop10pEdges_SpectralEnrichedGOterms.csv"

        # do enrichment for mcl clusters
        get_go_enrichments(graphfile, mcl_clustfile, mcl_outfile)
