from gothic import *

# - fdr max is 0.5? should be 1 [DONE]
# was doing fdr = numshufsignifs/numrealsignifs+numshufsignifs when should just be numshufsignifs/numrealsignifs
# cases = ['1a', '1b', '1c']
# for case in cases:
#     tpd_realvshuf_file = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_SignificantClustersRealVsShuffled.csv"
#     tpd_og_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_resultsDF.csv"
#     get_real_v_null_significants(tpd_og_results_file, startstopstep=(1, 0, -0.00001),
#                                  outfile=tpd_realvshuf_file, monotonic=True)

# - for fdr tpd figs, do with log y scale [DONE]
# did this in getBioinfoPaperFigs.py rather than here

# - look at 1c graph report for weirdness [ ]
cases = ['1a', '1b', '1c']
for case in cases:
    graphfile = "MEC" + case + "TADgraphPt0175cut100kwindow_oneminusw.gt"
    reportfile = "MEC" + case + "TADgraphPt0175cut100kwindow_GraphReport.txt"
    generate_graph_report(graphfile, outfile=reportfile, use_lcc=True)


# - heatmaps make them mirrored colours around 0.5 [ ]


# - redo plots so scales match [ ]


# - don't do linear analysis plot for mcl/spectral (doent make sense) [ ]


# - instead do comparison of linear distance distributions plots across all 3 methods [ ]


# - redo mcl/spectral heatmaps with pval threshold, and also plot all 1400 intersected terms [ ]


# - for spectral, why are there missing values in the heatmaps? [ ]


# ^ check when doing fullheatmaps [ ]


# - do cytoscape visualization of overlapping significant GOterms [ ]
#  get node lists for visualization
