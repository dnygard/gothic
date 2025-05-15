from gothic import *

for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
    print(case)
    run = case + "FixedBins500kb"
    infile = run + "_annotated.graphml"
    normfile = run + "_ICEnormed.gt"
    modfile = run + "_ICEnormedInverted.gt"
    gcoo = run + "_COO.tsv"
    normcoo = "_normedCOO.tsv"

    # balance and modify weights to by 1/normed contact frequency
    make_graph_coo(infile, outfile=gcoo)
    ice_balance(gcoo, normcoo)
    update_weights_from_coo(infile, normcoo, normfile)
    simple_invert_weights(normfile, wannotation="weight", outfile=modfile)

    # plot unfiltered edge weights
    regfgplot = run + "_Unfiltered_EdgeWeightsPlot.png"
    logfgplot = run + "_Unfiltered_EdgeWeightsLogPlot.png"
    plot_edge_weight_distribution(modfile, binsize=100, outfile=regfgplot, log=False)
    try:
        plot_edge_weight_distribution(modfile, binsize=100, outfile=logfgplot, log=True)
    except ValueError as ve:
        print("Can't make Logged plot: ")
        print(ve)

    #filter
    for flevel in [0.05, 0.1, 0.2, 0.5]:
        fpstring = "Top" + str(int(flevel*100)) + "p"
        fgfile = run + "_edgeModifiedFiltered" + fpstring + ".gt"
        regfgplot = run + "_" + fpstring + "_EdgeWeightsPlot.png"
        logfgplot = run + "_" + fpstring + "_EdgeWeightsLogPlot.png"
        print(fgfile)

        create_top_edges_graph(modfile, threshold=flevel, outfile=fgfile)

        plot_edge_weight_distribution(fgfile, binsize=100, outfile=regfgplot, log=False)
        try:
            plot_edge_weight_distribution(fgfile, binsize=100, outfile=logfgplot, log=True)
        except ValueError as ve:
            print("Can't make Logged plot: ")
            print(ve)
