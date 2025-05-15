from gothic import *

for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
    print(case)
    gfile = case + "FixedBins500kb_ICEnormed.gt"
    filteredfile = case + "FixedBins500kb_ICEnormedTop5p.gt"
    spectral_outfile = case + "_SpectralClusters.csv"
    MCL_outfile = case + "_MCLclusters.csv"
    spectral_enrichfile = case + "_SpectralEnrichments.csv"
    MCL_enrichfile = case + "_MCLenrichments.csv"

    # # filter
    # print("filtering...")
    # create_top_edges_graph(gfile, threshold=0.2, outfile=filteredfile, weights_are_distances=False)

    # cluster
    print("MCL clustering...")
    do_markov_clustering(filteredfile, MCL_outfile)
    print("Spectral clustering...")
    do_spectral_clustering(filteredfile, spectral_outfile)

    # get enrichments
    print("getting enriched terms...")
    get_go_enrichments(filteredfile, MCL_outfile, MCL_enrichfile)
    get_go_enrichments(filteredfile, spectral_outfile, spectral_enrichfile)
