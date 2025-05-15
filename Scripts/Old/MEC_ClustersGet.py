from gothic import *

# perform Spectral and MCL clustering on MEC replicates

for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
    print(case)
    gfile = case + "TADgraphPt0175cut100kwindow.gt"
    # filteredfile = case + "FixedBins500kb_ICEnormed_Interchromosomal.gt"
    spectral_outfile = case + "_SpectralClusters.csv"
    MCL_outfile = case + "_MCLclusters.csv"
    spectral_enrichfile = case + "_SpectralEnrichments.csv"
    MCL_enrichfile = case + "_MCLenrichments.csv"

    # # filter
    # print("filtering...")
    # extract_interchromosomal_subgraph(gfile, filteredfile, verbose=True)

    # # cluster
    # print("clustering...")
    # do_markov_clustering(filteredfile, MCL_outfile)
    # do_spectral_clustering(filteredfile, spectral_outfile)
    #
    # # get enrichments
    # get_go_enrichments(filteredfile, MCL_outfile, MCL_enrichfile)
    # get_go_enrichments(filteredfile, spectral_outfile, spectral_enrichfile)

    # cluster
    print("clustering...")
    do_markov_clustering(gfile, MCL_outfile)
    do_spectral_clustering(gfile, spectral_outfile)

    # get enrichments
    get_go_enrichments(gfile, MCL_outfile, MCL_enrichfile)
    get_go_enrichments(gfile, spectral_outfile, spectral_enrichfile)


