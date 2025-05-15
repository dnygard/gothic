from gothic import *

cases = ["1c"]
for case in cases:
    runname = "MEC" + case + "TADgraphPt0175cut100kwindow"

    # get real and shuffled TPDs and p-values for all go terms
    print("Getting real and shuffled TPDs and p-values...")
    graphname = runname + ".gt"
    newgraphname = runname + "_goshuffled.gt"
    matfilename = runname + "_distmatrix.tsv"
    gd = make_go_dict(graphname)  # need to make godict before doing all_go_tpd()
    shufgd = make_go_dict(newgraphname)  # need to make godict before doing all_go_tpd()
    all_go_tpd(gd, matfilename, runname + "_realTPDs.csv")
    all_go_tpd(shufgd, matfilename, runname + "_shuffledTPDs.csv")

    # get final dataframe with "nnodes", "tpd", "pval", "shuftpd", "shufpval" for all go terms
    print("Collecting final results...")
    graphname = runname + ".gt"
    distrnfile = runname + "_distributions.csv"
    resultscsvname = runname + "_resultsDF.csv"
    topgolistname = runname + "_topGOterms.csv"
    fdrthresholdsname = runname + "_FDRthresholds.csv"
    fdrplotname = runname + "_FDRpvalComparePlot.png"
    pvaldistrnplotname = runname + "_pvalHistogram.png"

    df = get_go_tpd_pvals(graphname, runname + "_realTPDs.csv", runname + "_shuffledTPDs.csv", distrnfile)
    df.to_csv(resultscsvname)

    # also get list of top go terms with FDR
    get_top_goterms(resultscsvname, outfile=topgolistname)
    get_real_v_null_significants(resultscsvname, outfile=fdrthresholdsname)

    # plots
    plot_real_shuffled_pval_distributions(resultscsvname, pvaldistrnplotname)
    plot_fdr_pval_histogram(resultscsvname, fdrplotname)

    # update results file with gonames and AB compartment statuses
    # take a list of GO ids, add GO terms to table
    obofile = "go-basic.obo"
    resultsfile = topgolistname
    outfile = runname + "_topGOtermsPLUSNAMESPLUSAB.csv"
    graphfile = runname + ".gt"

    # make goid:goname dictionary
    with open(obofile, "r") as gof:
        termdict = {}  # dictionary for mapping go ids to terms (key=ID, value=term)
        ids = []  # list to append ids and alt ids to
        for line in gof:
            splitline = line.split(": ")
            if splitline[0] == "name":
                goname = splitline[1].strip()
            if splitline[0] == "id":
                ids.append(splitline[1].strip())
            if splitline[0] == "altid":
                ids.append(splitline[1].strip())
            if splitline[0] == "def":
                for goid in ids:
                    termdict[goid] = goname
                ids = []  # empty list of ids

    for key in termdict.keys():
        print(key + "->" + termdict[key])
    print(len(termdict.keys()))
    print(termdict.keys())

    with open(resultsfile, "r") as rf:
        with open(outfile, "w") as outf:
            outf.write("goterm, pval, FDRatThisCutoff, goname, percentCompartmentA\n")
            for line in rf:
                if line.startswith(","):
                    continue
                splitline = line.split(",")
                goterm = splitline[1].strip()
                pval = splitline[3]
                fdr = splitline[4].strip()
                if goterm in termdict.keys():
                    goname = termdict[goterm]
                else:
                    goname = "Unknown"

                # get A/B compartment ratio of all nodes in cluster
                vlist = get_nodes_with_term(graphfile, goterm)
                g = load_graph(graphfile)
                ablist = [g.vp.ab_compartment[i] for i in vlist]
                abratio = sum(i == "A" for i in ablist) / len(ablist)

                outf.write(goterm + ", " + pval + ", " + fdr + ", " + goname + ", " + str(abratio) + "\n")


    # # create inverted weights graph for optics and mcl clustering
    #
    # # do spectral optics and mcl clustering
    # graphname = newgraphname
    # spectralresultsfile = runname + "_SpectralClusters.csv"
    # mclresultsfile = runname + "_MCLclusters.csv"
    #
    # # do clustering
    # do_markov_clustering(graphname, mclresultsfile, inflation=2.1)
    # do_spectral_clustering(graphname, spectralresultsfile)
    #
    # # do clustering enrichment
    # mcl_plotfile = runname + "_MCLenrichments.csv"
    # spectral_plotfile = runname + "_Spectralenrichments.csv"
    #
    # get_go_enrichments(graphname, mclresultsfile, mcl_plotfile)
    # get_go_enrichments(graphname, spectralresultsfile, spectral_plotfile)
