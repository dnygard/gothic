from gothic import *

cases = ["1a", "1b", "1c"]
for case in cases:
    runname = "MEC" + case + "TADgraphPt0175cut100kwindow"
    # do spectral optics and mcl clustering
    graphname = runname + ".gt"
    modgraph = runname + "_oneminusw.gt"
    spectralresultsfile = runname + "_SpectralClusters.csv"
    mclresultsfile = runname + "_MCLclusters.csv"


    # # create inverted weights graph for optics and mcl clustering
    # g = load_graph(graphname)
    # w_ep = g.new_edge_property("double")
    # for e in g.edges():
    #     w_ep[e] = 1 - g.ep.weight[e]
    #
    # g.ep.weight = w_ep
    #
    # g.save(modgraph)

    # do clustering
    do_markov_clustering(modgraph, mclresultsfile, inflation=1.3)
    do_spectral_clustering(modgraph, spectralresultsfile)

    # do clustering enrichment
    mcl_plotfile = runname + "_MCLenrichments.csv"
    spectral_plotfile = runname + "_Spectralenrichments.csv"

    get_go_enrichments(modgraph, mclresultsfile, mcl_plotfile)
    get_go_enrichments(modgraph, spectralresultsfile, spectral_plotfile)
