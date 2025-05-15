from gothic import *
from matplotlib.colors import LinearSegmentedColormap

cases = ["1a", "1b", "1c"]
plot2e = "MECTADgraphPt0175cut100kwindow_fig2e.png"
plot3d = "MECTADgraphPt0175cut100kwindow_fig3d.png"
plot4d = "MECTADgraphPt0175cut100kwindow_fig4d.png"

#constants
fdrthresh = 0.25
pvalthresh = 0.05
plotxranges = [0.1, 1]

# initialize dicts of dicts for all methods and cases
goabreplsdicttpd = {}  # dict of {"goname":{"1a":ABval}, {"1b":ABval}, {"1c":ABval}} fopr ab heatmap
for i in cases:
    goabreplsdicttpd[i] = {}  # make as dict of dicts (one dict for each case)

goabreplsdictmcl = {}  # dict of {"goname":{"1a":ABval}, {"1b":ABval}, {"1c":ABval}} fopr ab heatmap
for i in cases:
    goabreplsdictmcl[i] = {}  # make as dict of dicts (one dict for each case)

goabreplsdictspectral = {}  # dict of {"goname":{"1a":ABval}, {"1b":ABval}, {"1c":ABval}} fopr ab heatmap
for i in cases:
    goabreplsdictspectral[i] = {}  # make as dict of dicts (one dict for each case)

# initialize lists of lists for venn diagrams for all methods
vennlistlist = []
vennlistlist_mcl = []
vennlistlist_spectral = []

# define custom colormap for heatmaps
colors = ["#00085e", "#ffffff", "#00085e"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

for case in cases:
    # initialize all file names
    graphfile = "MEC" + case + "TADgraphPt0175cut100kwindow_oneminusw.gt"
    #tpd_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_topGOterms.csv"
    tpd_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_SupplementalTable.csv"
    #tpd_og_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_ShufvRealFDRtable.csv"
    tpd_revigofile = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_RevigoList.txt"
    tpd_og_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_resultsDF.csv"
    tpd_realvshuf_file = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_SignificantClustersRealVsShuffled.csv"
    tpd_vennlist = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_vennlist.txt"
    mcl_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_SupplementalTable.csv"
    mcl_enrichment_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLenrichments.csv"
    mcl_cluster_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters.csv"
    mcl_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters_Min3Filtered.csv"
    mcl_vennlist = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_vennlist.txt"
    mcl_plotfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_LinDistPlot.png"
    mcl_revigofile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_RevigoList.txt"
    spectral_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_SupplementalTable.csv"
    spectral_enrichment_file = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectralenrichments.csv"
    spectral_cluster_file = "MEC" + case + "TADgraphPt0175cut100kwindow_SpectralClusters.csv"
    spectral_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_SpectralClusters_Min3Filtered.csv"
    spectral_vennlist = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_vennlist.txt"
    spectral_plotfile = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_LinDistPlot.png"
    spectral_revigofile = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_RevigoList.txt"
    gaf_file = "goa_human.gaf"
    gencode_file = "gencode_pcgenes.csv"
    golist_file = "MEC" + case + "TADgraphPt0175cut100kwindow_AnnotatedGOterms.txt"
    goobo_file = "go-basic.obo"
    plot2a = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2a.png"
    plot2b = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2b.png"
    plot2c = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2c.png"
    plot2d = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2d.png"
    plot3a = "MEC" + case + "TADgraphPt0175cut100kwindow_fig3a.png"
    plot3b = "MEC" + case + "TADgraphPt0175cut100kwindow_fig3b.png"
    plot3c = "MEC" + case + "TADgraphPt0175cut100kwindow_fig3c.png"
    plot3e = "MEC" + case + "TADgraphPt0175cut100kwindow_fig3e.png"
    plot4a = "MEC" + case + "TADgraphPt0175cut100kwindow_fig4a.png"
    plot4b = "MEC" + case + "TADgraphPt0175cut100kwindow_fig4b.png"
    plot4c = "MEC" + case + "TADgraphPt0175cut100kwindow_fig4c.png"
    plot4e = "MEC" + case + "TADgraphPt0175cut100kwindow_fig4e.png"

    #### Fig 1 pipeline/workflow done in photoshop/illustrator

    #### Fig 2: TPD for all 3 Mammary Epithelial cell replicates
        # Figure 2: A) B) C) D) E)

    # read in results file as lists
    with open(tpd_results_file, "r") as rf:
        gotermlist = []
        gonamelist = []
        nodeswithtermlist = []
        chromswithtermlist = []
        tpdlist = []
        pvallist = []
        fdrlist = []
        abratiolist = []
        avglineardistancelist = []
        shuftpdlist = []
        shufpvallist = []

        for line in rf:
            if not line.startswith("GOterm"):  # skip header line
                try:
                    splitline = line.strip().split(",")
                    gotermlist.append(splitline[0])
                    nodeswithtermlist.append(splitline[1])
                    chromswithtermlist.append(splitline[2])
                    tpdlist.append(float(splitline[3]))
                    pvallist.append(float(splitline[4]))
                    fdrlist.append(float(splitline[5]))
                    abratiolist.append(float(splitline[6]))
                    avglineardistancelist.append(float(splitline[7]))
                    shuftpdlist.append(float(splitline[8]))
                    shufpvallist.append(float(splitline[9]))
                    gonamelist.append(splitline[10])
                except IndexError:
                    print("INDEXERROR: " + line)

    for xrange in plotxranges:  # redo figures for different x limits
        plot2a = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2a.png"
        plot2b = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2b.png"
        plot2c = "MEC" + case + "TADgraphPt0175cut100kwindow_fig2c.png"
        plot2a = plot2a.split(".")[0] + str(xrange).replace(".", "pt") + ".png"
        plot2b = plot2b.split(".")[0] + str(xrange).replace(".", "pt") + ".png"
        plot2c = plot2c.split(".")[0] + str(xrange).replace(".", "pt") + ".png"

        # A) Scatter plot of #of GO terms (y-axis) against FDR (x-axis) [ ]
        get_real_v_null_significants(tpd_og_results_file, startstopstep=(xrange, 0, -xrange/10000),
                                     outfile=tpd_realvshuf_file, monotonic=True)
        # get ymax for plotting
        df = pd.read_csv(tpd_realvshuf_file, sep=", ")  # read input csv into data frame
        df = df.astype({"pvalThreshold": "float", "numReal": "int", "numShuf": "int", "fdr": "float"})
        numtermslist = []
        for index, row in df.iterrows():
            numtermslist.append(row["numReal"])
        numtermslist = [int(x) for x in numtermslist]
        yrange = max(numtermslist)

        plot_passing_terms_v_fdr_fig(tpd_realvshuf_file, outfile=plot2a, xlim=xrange, ylim=yrange)

        # redo plots with log yscale
        logplot2a = plot2a.split(".")[0] + "_LogYscale.png"
        plot_passing_terms_v_fdr_fig(tpd_realvshuf_file, outfile=logplot2a, xlim=xrange, ylim=yrange, logy=True)

        # B) FDR (y-axis) against p-values (x-axis) [ ]
        #plot_fdr_pval_histogram(results_file, outfile, stepsize=0.001, log_axis=True)
        stepsize = xrange/100
        logaxis = True
        pylist = [float(x) for x in pvallist]  # convert to numeric
        shufpylist = [float(x) for x in shufpvallist]  # convert to numeric
        xlist = list(np.arange(stepsize, xrange + stepsize, stepsize))

        if logaxis:
            #logbins = np.logspace(np.log10(xlist[0]), np.log10(xlist[-1]), len(xlist))
            logbins = xlist
        else:
            logbins = xlist

        plt.clf()
        # Create histogram
        plt.hist(pylist, bins=logbins, alpha=0.5, label='real graph', color='blue', cumulative=True)
        plt.hist(shufpylist, bins=logbins, alpha=0.5, label='null graph', color='black', cumulative=True)

        # Add labels and title
        plt.xlabel('p-value threshold', fontsize=20)
        plt.ylabel('# of GO terms with\n p-value < threshold', fontsize=20)
        plt.xlim(xlist[0], xlist[-1])  # why does setting this result in such a strange plot?
        if logaxis:
            plt.yscale('log')
            #plt.xscale('log')  # TODO remove this line
        plt.legend(loc='upper left')
        # Save the figure
        plt.savefig(plot2b, bbox_inches='tight')
        plt.close()

        # C) Distribution of p-values for randomized annotation and real ones [ ]
        #plot_real_shuffled_pval_distributions(results_file, outfile, stepsize=0.001)
        stepsize = xrange/100
        pylist = [float(x) for x in pvallist]  # convert to numeric
        shufpylist = [float(x) for x in shufpvallist]  # convert to numeric
        xlist = list(np.arange(stepsize, xrange + stepsize, stepsize))

        plt.clf()
        # Create histogram
        plt.hist(pylist, bins=xlist, alpha=0.5, label='real graph', color='blue', cumulative=False)
        plt.hist(shufpylist, bins=xlist, alpha=0.5, label='null graph', color='black', cumulative=False)

        # Add labels and title
        plt.xlabel('p-value', fontsize=20)
        plt.ylabel('# of GO terms in\n p-value range', fontsize=20)
        plt.xlim(xlist[0], xlist[-1])
        plt.legend(loc='upper right')
        # Save the figure
        plt.savefig(plot2c, bbox_inches='tight')
        plt.close()

    # D) Scatterplot of linear distance (y-axis) against ranking of GO terms by TPD significance (p-values) [ ]
    #get_linear_distances_tpd(resultsfile, gafname, gencodefile, outfile=None, nsamples=1000000, verbose=True, plot=True)
    # replace -1 distances with None
    ilist = []
    for i in range(len(avglineardistancelist)):
        if avglineardistancelist[i] != -1:
            ilist.append(i)  # list of indices with actual linear distance values
    reallindistlist = [avglineardistancelist[x] for x in ilist]

    # get mean and standard deviation for random gene pairs
    mean = sum(reallindistlist) / len(reallindistlist)
    variance = sum([((x - mean) ** 2) for x in reallindistlist]) / len(reallindistlist)
    stdev = variance ** 0.5
    # Create a scatter plot
    plt.figure().set_figwidth(15)
    plt.scatter(ilist, reallindistlist, s=0.5, label='Linear Distance vs GO Ranking', color='black')

    # Set labels for the x and y axes
    plt.xlabel('Ranking of GO term by TPD significance')
    plt.ylabel('Linear distance (bp)')

    corr, _ = pearsonr(ilist, reallindistlist)
    plt.title("R^2 = " + str(corr), loc='left')
    # add lines for avg linear distance and standard deviations
    meanlabel = "mean = " + str(mean)
    upperlabel = "stddev = " + str(stdev)
    plt.axhline(y=mean, color='red', linestyle='-', label=meanlabel)  # mean
    plt.axhline(y=mean + stdev, color='red', linestyle='--', label=upperlabel)  # stddev upper bound
    plt.axhline(y=mean - stdev, color='red', linestyle='--')  # stddev lower bound
    plt.legend()
    plt.savefig(plot2d)

    # alternative way to do linear comparison where we just compare the linear distance distributions between clusters and random
    print("TPD " + case)
    regions_file = "MEC" + case + "TADgraphPt0175cut100kwindow.tads"
    tpd_lindisthist = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_LinearDistanceDistributionsPlot.png"
    g = load_graph(graphfile)
    # create node position dict of dicts from TAD regions file
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    posdict = {}
    for chrom in chrlist:  # initialize dict of lists so posdict[chr1][0] is tuple of positions for first node of that chrom
        posdict[chrom] = []
    with open(regions_file, "r") as f:
        for line in f:
            chrom, start, stop = line.split()
            posdict[chrom].append((start, stop))  # pos is stored as a (start, stop) tuple

    # get list of clusters
    clustlist = []
    print("getting clusters...")
    for term in tqdm(gotermlist[0:99], total=100):
        clustlist.append(get_nodes_with_term(g, term))
    # for each cluster, get all pairwise distances for nodes in cluster
    distlist = []
    print("getting distances...")
    for clust in tqdm(clustlist, total=len(clustlist)):
        vcombs = list(combinations(clust, 2))
        for n1, n2 in vcombs:
            vchrom1, vnode1 = g.vp.vname[n1].split(":")
            vchrom2, vnode2 = g.vp.vname[n2].split(":")
            vnode1 = int(vnode1)
            vnode2 = int(vnode2)
            start1, stop1 = posdict[vchrom1][vnode1]
            start2, stop2 = posdict[vchrom2][vnode2]
            start1 = int(start1)
            start2 = int(start2)
            stop1 = int(stop1)
            stop2 = int(stop2)

            if vchrom1 == vchrom2:
                if start1 > stop2 and stop1 > stop2:  # if 1 is later in chrom than 2
                    dist = start1 - stop2
                elif start1 < start2 and stop1 < start2:  # if 2 is later in chrom than 1
                    dist = start2 - stop1
                else:
                    dist = 0  # if any overlap dist is 0
                distlist.append(dist)

    # get dists for random node pairs belonging to the same (random) chromosome to compare against
    randdistlist = []
    print("getting distances for random...")
    for i in tqdm(range(len(distlist)), total=len(distlist)):
        # get random chromosome
        rchrom = random.choice(chrlist)
        while len(posdict[rchrom]) < 2:  # ensure we are picking from a list with at least two entries
            rchrom = random.choice(chrlist)
        start1, stop1 = random.choice(posdict[rchrom])
        start2, stop2 = random.choice(posdict[rchrom])
        while start1 == start2:
            start2, stop2 = random.choice(posdict[rchrom])  # ensure we pick two different nodes
        start1 = int(start1)
        start2 = int(start2)
        stop1 = int(stop1)
        stop2 = int(stop2)

        if start1 > stop2 and stop1 > stop2:  # if 1 is later in chrom than 2
            dist = start1 - stop2
        elif start1 < start2 and stop1 < start2:  # if 2 is later in chrom than 1
            dist = start2 - stop1
        else:
            dist = 0  # if any overlap dist is 0
        randdistlist.append(dist)

    if distlist and randdistlist:
        print("plotting...")
        # plot distributions as overlapping bar plots
        xmax = max(max(distlist), max(randdistlist))  # get max value for x axis
        xstep = int(xmax / 100)  # bin size to make sure there are always 100 bins
        xlist = np.arange(0, int(xmax), xstep)  # get list of bins for plot
        # plot
        plt.clf()
        # Create histogram
        plt.hist(distlist, bins=xlist, alpha=0.5, label='cluster pairwise linear distances', color='blue')
        plt.hist(randdistlist, bins=xlist, alpha=0.5, label='random pairwise linear distances', color='black')
        # Add labels and title
        plt.xlabel('linear distance (bp)', fontsize=20)
        plt.ylabel('number of occurences', fontsize=20)
        # plt.xlim(xlist[0], xlist[-1])  # why does setting this result in such a strange plot?
        # if log_axis:
        #     plt.yscale('log')
        #     plt.xscale('log')  # TODO remove this line
        plt.legend(loc='upper left')
        # Save the figure
        plt.savefig(tpd_lindisthist, bbox_inches='tight')
        plt.close()


    # E) A/B compartments analysis with the ratio in a heatmap where all 3 replicates are combined (next to each other as columns) on the heatmap [ ]
    # Store output in abtpdlistlist to be plotted after looping through all replicates
    # need to store both GO names and associated ab ratios (tuple?)
    for i in range(len(abratiolist)):
        goabreplsdicttpd[case][gotermlist[i]] = abratiolist[i]  # store ab values in dict of dicts

    # F) Venn Diagram of significant GO terms (need to choose an FDR based on above results) where each replicate is a circle. [ ]
    # Supplemental 1: Revigo figures of significant GO terms (need to choose an FDR based on above results) [ ]
        # DONE VIA GUI, NOT DONE PROGRAMATICALLY HERE
    # build list of lists to be sent to venn diagram software
    vennlist = []
    with open(tpd_vennlist, "w") as f:
        for i in range(len(fdrlist)):
            if fdrlist[i] <= fdrthresh:
                vennlist.append(gotermlist[i])
                f.write(gotermlist[i] + "\n")
    vennlistlist.append(vennlist)
#
#     # 2supp) Revigo figures of significant GO terms (need to choose an FDR-adjusted enrichment p-value based on above results) [ ]
    # write lists of go terms to be sent to revigo
    with open(tpd_revigofile, "w") as f:
        for i in range(len(fdrlist)):
            if pvallist[i] <= pvalthresh:
                if fdrlist[i] <= fdrthresh:
                    f.write(gotermlist[i] + " " + str(pvallist[i]) + "\n")

    # clear lists for tpd results from memory
    del gotermlist
    del gonamelist
    del nodeswithtermlist
    del chromswithtermlist
    del tpdlist
    del pvallist
    del fdrlist
    del abratiolist
    del avglineardistancelist


    #### Fig 3: MCL and GO analysis for all 3 Mammary Epithelial cell replicates [ ]
    # read in results file as lists
    with open(mcl_results_file, "r") as rf:
        gotermlist = []
        gonamelist = []
        nodeswithterminclusterlist = []
        numberofnodeswithtermingraphlist = []
        numberofnodeswithterminclusterlist = []
        chromswithtermlist = []
        chromswithterminclusterlist = []
        pvallist = []
        abratiolist = []
        avglineardistancelist = []

        for line in rf:
            if not line.startswith("GOterm"):  # skip header line
                splitline = line.strip().split(",")
                gotermlist.append(splitline[0])
                nodeswithterminclusterlist.append(splitline[1])
                numberofnodeswithtermingraphlist.append(splitline[2])
                numberofnodeswithterminclusterlist.append(splitline[3])
                chromswithtermlist.append(splitline[4])
                chromswithterminclusterlist.append(splitline[5])
                pvallist.append(float(splitline[6]))
                abratiolist.append(float(splitline[7]))
                avglineardistancelist.append(float(splitline[8]))
                gonamelist.append(splitline[9])

    # 3a) Distribution of cluster size (showing how many cluster with size 1, 2, 3, ……) [ ]
    # count cluster sizes
    sizelist = []
    with open(mcl_cluster_file, "r") as cf:
        for line in cf:
            sizelist.append(len(line.split()))

    pylist = sizelist

    plt.clf()
    # Create histogram
    plt.hist(pylist, bins=max(pylist), alpha=0.5, label='real graph', color='black', cumulative=False)

    # Add labels and title
    plt.xlabel('number of nodes in cluster', fontsize=20)
    plt.ylabel('count', fontsize=20)
    # Save the figure
    plt.savefig(plot3a, bbox_inches='tight')
    plt.close()

    # 3b) Distribution of the average linear distances of the set of nodes in a cluster annotated with significant GO terms.
    # Each count in the histogram will correspond to 1 average of the linear distances of all pairs of nodes in a cluster
    # that are annotated with a significant GO term, with a line indicating the average linear distance in the network.
    # We ignore pairs of nodes in different chromosomes when calculating the average. [✓]
    mcl_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLenrichments.csv"
    mcl_cluster_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters.csv"
    mcl_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters_Min3Filtered.csv"
    mcl_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_SupplementalTable.csv"
    mcl_plotfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_LinDistPlot.png"
    do_linear_analysis_mcl(graphfile, mcl_enrichment_file, gencode_file, plotfile=plot3b)
    # ^deprecated because it makes no sense to do linear analysis this way for MCL/Spectral
    # instead, just do distribution plots (histograms) of pairwise linear distances (also do for TPD)
    regions_file = "MEC" + case + "TADgraphPt0175cut100kwindow.tads"
    mcl_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters.csv"
    mcl_lindisthist = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_LinearDistanceDistributionsPlot.png"
    g = load_graph(graphfile)
    # create node position dict of dicts from TAD regions file
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    posdict = {}
    for chrom in chrlist:  # initialize dict of lists so posdict[chr1][0] is tuple of positions for first node of that chrom
        posdict[chrom] = []
    with open(regions_file, "r") as f:
        for line in f:
            chrom, start, stop = line.split()
            posdict[chrom].append((start, stop))  # pos is stored as a (start, stop) tuple

    # get list of clusters
    clustlist = []
    with open(mcl_cluster_outfile, "r") as f:
        for line in f:
            checklist = [int(x.rstrip(",")) for x in line.strip().split(", ")]
            if len(checklist) >= 3:
                clustlist.append(checklist)  # convert line to list of ints
    # for each cluster, get all pairwise distances for nodes in cluster
    distlist = []
    for clust in clustlist:
        vcombs = list(combinations(clust, 2))
        for n1, n2 in vcombs:
            vchrom1, vnode1 = g.vp.vname[n1].split(":")
            vchrom2, vnode2 = g.vp.vname[n2].split(":")
            vnode1 = int(vnode1)
            vnode2 = int(vnode2)
            start1, stop1 = posdict[vchrom1][vnode1]
            start2, stop2 = posdict[vchrom2][vnode2]
            start1 = int(start1)
            start2 = int(start2)
            stop1 = int(stop1)
            stop2 = int(stop2)

            if vchrom1 == vchrom2:
                if start1 > stop2 and stop1 > stop2:  # if 1 is later in chrom than 2
                    dist = start1 - stop2
                elif start1 < start2 and stop1 < start2:  # if 2 is later in chrom than 1
                    dist = start2 - stop1
                else:
                    dist = 0  # if any overlap dist is 0
                distlist.append(dist)

    # get dists for random node pairs belonging to the same (random) chromosome to compare against
    randdistlist = []
    for i in range(len(distlist)):
        # get random chromosome
        rchrom = random.choice(chrlist)
        while len(posdict[rchrom]) < 2:  # ensure we are picking from a list with at least two entries
            rchrom = random.choice(chrlist)
        start1, stop1 = random.choice(posdict[rchrom])
        start2, stop2 = random.choice(posdict[rchrom])
        while start1 == start2:
            start2, stop2 = random.choice(posdict[rchrom])  # ensure we pick two different nodes
        start1 = int(start1)
        start2 = int(start2)
        stop1 = int(stop1)
        stop2 = int(stop2)

        if start1 > stop2 and stop1 > stop2:  # if 1 is later in chrom than 2
            dist = start1 - stop2
        elif start1 < start2 and stop1 < start2:  # if 2 is later in chrom than 1
            dist = start2 - stop1
        else:
            dist = 0  # if any overlap dist is 0
        randdistlist.append(dist)

    if distlist and randdistlist:
        # plot distributions as overlapping bar plots
        xmax = max(max(distlist), max(randdistlist))  # get max value for x axis
        xstep = int(xmax/100)  # bin size to make sure there are always 100 bins
        xlist = np.arange(0, int(xmax), xstep)  # get list of bins for plot
        # plot
        plt.clf()
        # Create histogram
        plt.hist(distlist, bins=xlist, alpha=0.5, label='cluster pairwise linear distances', color='blue')
        plt.hist(randdistlist, bins=xlist, alpha=0.5, label='random pairwise linear distances', color='black')
        # Add labels and title
        plt.xlabel('linear distance (bp)', fontsize=20)
        plt.ylabel('number of occurences', fontsize=20)
        # plt.xlim(xlist[0], xlist[-1])  # why does setting this result in such a strange plot?
        # if log_axis:
        #     plt.yscale('log')
        #     plt.xscale('log')  # TODO remove this line
        plt.legend(loc='upper left')
        # Save the figure
        plt.savefig(mcl_lindisthist, bbox_inches='tight')
        plt.close()

    # 3c) Revigo figures of significant GO terms (need to choose an FDR-adjusted enrichment p-value based on above results) [ ]
    # write lists of go terms to be sent to revigo
    with open(mcl_revigofile, "w") as f:
        for i in range(len(pvallist)):
            if pvallist[i] <= pvalthresh:
                f.write(gotermlist[i] + " " + str(pvallist[i]) + "\n")

    # 3d) A/B compartments analysis with the ratio in a heatmap where all 3 replicates are combined (next to each other as
    # columns) on the heatmap for significant GO terms. Ratio is calculated only on the set of nodes in a cluster annotated with the significant GO terms. [ ]
    # DONE LATER SO ALL 3 REPLICATES CAN BE COMBINED
    # Store output in abtpdlistlist to be plotted after looping through all replicates
    # need to store both GO names and associated ab ratios (tuple?)
    for i in range(len(abratiolist)):
        goabreplsdictmcl[case][gotermlist[i]] = abratiolist[i]  # store ab values in dict of dicts

    # 3e) Venn Diagram of significant GO terms (need to choose an FDR-adjusted enrichment p-value based on above results)
    # where each replicate is a circle. [ ]

    # build list of lists to be sent to venn diagram software
    vennlist = []
    with open(mcl_vennlist, "w") as f:
        for i in range(len(pvallist)):
            if pvallist[i] <= pvalthresh:
                vennlist.append(gotermlist[i])
                f.write(gotermlist[i] + "\n")
    vennlistlist_mcl.append(vennlist)

    # Supplementary Files: Go analysis for each cluster of MCL [ ]
    # Supplementary Table: Where each line is a GO terms and columns feature, GO ID, GO name, GO size, Cluster ID,
    # Number of nodes in the cluster, Number of nodes with GO annotation in the network, Number of nodes with GO annotation
    # in the cluster, Enrichment p-value, FDR-adjusted enrichment p-value, linear distance of the nodes with the GO term in
    # the cluster, A/B ratio of the nodes with the GO term in the cluster. [ ]

    # clear mcl lists from memory
    del gotermlist
    del gonamelist
    del nodeswithterminclusterlist
    del numberofnodeswithtermingraphlist
    del numberofnodeswithterminclusterlist
    del chromswithtermlist
    del chromswithterminclusterlist
    del pvallist
    del abratiolist
    del avglineardistancelist


    # #### Fig 4: Spectral clustering and GO analysis for all 3 Mammary Epithelial cell replicates
    # read in results file as lists
    with open(spectral_results_file, "r") as rf:
        gotermlist = []
        gonamelist = []
        nodeswithterminclusterlist = []
        numberofnodeswithtermingraphlist = []
        numberofnodeswithterminclusterlist = []
        chromswithtermlist = []
        chromswithterminclusterlist = []
        pvallist = []
        abratiolist = []
        avglineardistancelist = []

        for line in rf:
            if not line.startswith("GOterm"):  # skip header line
                splitline = line.strip().split(",")
                gotermlist.append(splitline[0])
                nodeswithterminclusterlist.append(splitline[1])
                numberofnodeswithtermingraphlist.append(splitline[2])
                numberofnodeswithterminclusterlist.append(splitline[3])
                chromswithtermlist.append(splitline[4])
                chromswithterminclusterlist.append(splitline[5])
                pvallist.append(float(splitline[6]))
                abratiolist.append(float(splitline[7]))
                avglineardistancelist.append(float(splitline[8]))
                gonamelist.append(splitline[9])
#
    # 4a) Distribution of cluster size (showing how many cluster with size 1, 2, 3, ……) [ ]
    # count cluster sizes
    sizelist = []
    with open(spectral_cluster_file, "r") as cf:
        for line in cf:
            sizelist.append(len(line.split(",")))

    pylist = sizelist

    plt.clf()
    # Create histogram
    plt.hist(pylist, bins=max(pylist), alpha=0.5, label='real graph', color='black', cumulative=False)

    # Add labels and title
    plt.xlabel('number of nodes in cluster', fontsize=20)
    plt.ylabel('count', fontsize=20)
    # Save the figure
    plt.savefig(plot4a, bbox_inches='tight')
    plt.close()

    # 4b) Distribution of the average linear distances of the set of nodes in a cluster annotated with significant GO terms.
    # Each count in the histogram will correspond to 1 average of the linear distances of all pairs of nodes in a cluster
    # that are annotated with a significant GO term, with a line indicating the average linear distance in the network.
    # We ignore pairs of nodes in different chromosomes when calculating the average. [✓]
    # mcl_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLenrichments.csv"
    # mcl_cluster_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters.csv"
    # mcl_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters_Min3Filtered.csv"
    # mcl_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_SupplementalTable.csv"
    # mcl_plotfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_LinDistPlot.png"
    do_linear_analysis_mcl(graphfile, spectral_enrichment_file, gencode_file, plotfile=plot4b)

    regions_file = "MEC" + case + "TADgraphPt0175cut100kwindow.tads"
    spectral_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_SpectralClusters.csv"
    spectral_lindisthist = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_LinearDistanceDistributionsPlot.png"
    g = load_graph(graphfile)
    # create node position dict of dicts from TAD regions file
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    posdict = {}
    for chrom in chrlist:  # initialize dict of lists so posdict[chr1][0] is tuple of positions for first node of that chrom
        posdict[chrom] = []
    with open(regions_file, "r") as f:
        for line in f:
            chrom, start, stop = line.split()
            posdict[chrom].append((start, stop))  # pos is stored as a (start, stop) tuple

    # get list of clusters
    clustlist = []
    with open(spectral_cluster_outfile, "r") as f:
        for line in f:
            checklist = [int(x) for x in line.strip().split(",")]
            if len(checklist) >= 3:
                clustlist.append(checklist)  # convert line to list of ints
    # for each cluster, get all pairwise distances for nodes in cluster
    distlist = []
    for clust in clustlist:
        vcombs = list(combinations(clust, 2))
        for n1, n2 in vcombs:
            vchrom1, vnode1 = g.vp.vname[n1].split(":")
            vchrom2, vnode2 = g.vp.vname[n2].split(":")
            vnode1 = int(vnode1)
            vnode2 = int(vnode2)
            start1, stop1 = posdict[vchrom1][vnode1]
            start2, stop2 = posdict[vchrom2][vnode2]
            start1 = int(start1)
            start2 = int(start2)
            stop1 = int(stop1)
            stop2 = int(stop2)

            if vchrom1 == vchrom2:
                if start1 > stop2 and stop1 > stop2:  # if 1 is later in chrom than 2
                    dist = start1 - stop2
                elif start1 < start2 and stop1 < start2:  # if 2 is later in chrom than 1
                    dist = start2 - stop1
                else:
                    dist = 0  # if any overlap dist is 0
                distlist.append(dist)

    # get dists for random node pairs belonging to the same (random) chromosome to compare against
    randdistlist = []
    for i in range(len(distlist)):
        # get random chromosome
        rchrom = random.choice(chrlist)
        while len(posdict[rchrom]) < 2:  # ensure we are picking from a list with at least two entries
            rchrom = random.choice(chrlist)
        start1, stop1 = random.choice(posdict[rchrom])
        start2, stop2 = random.choice(posdict[rchrom])
        while start1 == start2:
            start2, stop2 = random.choice(posdict[rchrom])  # ensure we pick two different nodes
        start1 = int(start1)
        start2 = int(start2)
        stop1 = int(stop1)
        stop2 = int(stop2)

        if start1 > stop2 and stop1 > stop2:  # if 1 is later in chrom than 2
            dist = start1 - stop2
        elif start1 < start2 and stop1 < start2:  # if 2 is later in chrom than 1
            dist = start2 - stop1
        else:
            dist = 0  # if any overlap dist is 0
        randdistlist.append(dist)

    if distlist and randdistlist:
        # plot distributions as overlapping bar plots
        xmax = max(max(distlist), max(randdistlist))  # get max value for x axis
        xstep = int(xmax / 100)  # bin size to make sure there are always 100 bins
        xlist = np.arange(0, int(xmax), xstep)  # get list of bins for plot
        # plot
        plt.clf()
        # Create histogram
        plt.hist(distlist, bins=xlist, alpha=0.5, label='cluster pairwise linear distances', color='blue')
        plt.hist(randdistlist, bins=xlist, alpha=0.5, label='random pairwise linear distances', color='black')
        # Add labels and title
        plt.xlabel('linear distance (bp)', fontsize=20)
        plt.ylabel('number of occurences', fontsize=20)
        # plt.xlim(xlist[0], xlist[-1])  # why does setting this result in such a strange plot?
        # if log_axis:
        #     plt.yscale('log')
        #     plt.xscale('log')  # TODO remove this line
        plt.legend(loc='upper left')
        # Save the figure
        plt.savefig(spectral_lindisthist, bbox_inches='tight')
        plt.close()

    # 4c) Revigo figures of significant GO terms (need to choose an FDR-adjusted enrichment p-value based on above results) [ ]
    # write lists of go terms to be sent to revigo
    with open(spectral_revigofile, "w") as f:
        for i in range(len(pvallist)):
            if pvallist[i] <= pvalthresh:
                f.write(gotermlist[i] + " " + str(pvallist[i]) + "\n")

    # 4d) A/B compartments analysis with the ratio in a heatmap where all 3 replicates are combined (next to each other as
    # columns) on the heatmap for significant GO terms. Ratio is calculated only on the set of nodes in a cluster annotated with the significant GO terms. [ ]
    # DONE LATER SO ALL 3 REPLICATES CAN BE COMBINED
    # Store output in abtpdlistlist to be plotted after looping through all replicates
    # need to store both GO names and associated ab ratios (tuple?)
    for i in range(len(abratiolist)):
        goabreplsdictspectral[case][gotermlist[i]] = abratiolist[i]  # store ab values in dict of dicts

    # 4e) Venn Diagram of significant GO terms (need to choose an FDR-adjusted enrichment p-value based on above results)
    # where each replicate is a circle. [ ]

    # build list of lists to be sent to venn diagram software
    vennlist = []
    with open(spectral_vennlist, "w") as f:
        for i in range(len(pvallist)):
            if pvallist[i] <= pvalthresh:
                vennlist.append(gotermlist[i])
                f.write(gotermlist[i] + "\n")
    vennlistlist_spectral.append(vennlist)

    # Supplementary Files: Go analysis for each cluster of MCL [ ]
    # Supplementary Table: Where each line is a GO terms and columns feature, GO ID, GO name, GO size, Cluster ID,
    # Number of nodes in the cluster, Number of nodes with GO annotation in the network, Number of nodes with GO annotation
    # in the cluster, Enrichment p-value, FDR-adjusted enrichment p-value, linear distance of the nodes with the GO term in
    # the cluster, A/B ratio of the nodes with the GO term in the cluster. [ ]

    # clear mcl lists from memory
    del gotermlist
    del gonamelist
    del nodeswithterminclusterlist
    del numberofnodeswithtermingraphlist
    del numberofnodeswithterminclusterlist
    del chromswithtermlist
    del chromswithterminclusterlist
    del pvallist
    del abratiolist
    del avglineardistancelist


#### Fig 5: Overlap between all 3 clustering methods

# Venn Diagram of intersections of the above 3 Venn Diagrams [ ]

Maybe upset plot of all GO identified with all replicates being independent. [ ]

Visualization of a few interesting GO terms perhaps some that are detected by all 3 methods [ ]

# plot 2E ab heatmaps
# get list of goterms for heatmap by taking intersections of vennlists
combs = combinations(vennlistlist, 2)
intxlist = []
for i, j in combs:
    intx = list(set(i) & set(j))
    intxlist.extend(intx)
intxset = list(set(intxlist))  # final list of goterms to be used

# create ndarray for storing ab vals
arrayshape = (len(intxset), len(cases))
abvalsarray = np.ndarray(arrayshape)
# populate array with abvals
for gotermi in range(len(intxset)):
    for casei in range(len(cases)):
        try:
            abvalsarray[gotermi, casei] = goabreplsdicttpd[cases[casei]][intxset[gotermi]]  # pulls ab val from dict
        except KeyError as ke:
            print(ke)
            abvalsarray[gotermi, casei] = np.nan  # if for some reason val not found, print error and fill with no val

print(abvalsarray.T)

# add GO names to term list for use as axis labels
# make goid:goname dictionary
gonames_list = []
print("creating goid:goname dictionary...")
with open(goobo_file, "r") as gof:
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
for go in intxset:  # just use first dict for iterating through so order is consistent
    goname = termdict[go]  # fetch name of goterm from obo
    if len(goname) > 37:
        goname = goname[:38] + "..."  # truncate goname if it is too long so it fits on plot
    gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later

# create heatmap
fig, ax = plt.subplots()
im = ax.imshow(abvalsarray, cmap=custom_cmap)

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(cases)), labels=cases)
ax.set_yticks(np.arange(len(gonames_list)), labels=gonames_list)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

cbar = plt.colorbar(im)
cbar.ax.set_ylabel("<- 100% A compartment           100% B compartment ->",
                   rotation=270, horizontalalignment="center", labelpad=20)
# # Loop over data dimensions and create text annotations.
# for i in range(len(gonames_list)):
#     for j in range(len(cases)):
#         text = ax.text(j, i, abvalsarray[i, j],
#                        ha="center", va="center", color="w")

fig.tight_layout()
plt.savefig(plot2e)
plt.clf()


# plot Fig3D ab heatmaps
# get list of goterms for heatmap by taking intersections of vennlists
combs = combinations(vennlistlist_mcl, 2)
intxlist = []
for i, j in combs:
    intx = list(set(i) & set(j))
    intxlist.extend(intx)
intxset = list(set(intxlist))  # final list of goterms to be used

print(len(intxset))
if len(intxset) > 30:
    intxset = intxset[0:29]
print(len(intxset))

# create ndarray for storing ab vals
arrayshape = (len(intxset), len(cases))
abvalsarray = np.ndarray(arrayshape)
# populate array with abvals
for gotermi in range(len(intxset)):
    for casei in range(len(cases)):
        try:
            print(cases[casei] + " - " + intxset[gotermi])
            abvalsarray[gotermi, casei] =         mcl[cases[casei]][intxset[gotermi]]  # pulls ab val from dict
        except KeyError as ke:
            print(ke)
            abvalsarray[gotermi, casei] = np.nan  # if for some reason val not found, print error and fill with no val

print(abvalsarray.T)

# add GO names to term list for use as axis labels
# make goid:goname dictionary
gonames_list = []
print("creating goid:goname dictionary...")
with open(goobo_file, "r") as gof:
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

for go in intxset:  # just use first dict for iterating through so order is consistent
    try:
        goname = termdict[go]  # fetch name of goterm from obo
        if len(goname) > 37:
            goname = goname[:38] + "..."  # truncate goname if it is too long so it fits on plot
        gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later
    except KeyError as ke:
        gonames_list.append(go + " - " + "Unknown")

# create heatmap
fig, ax = plt.subplots()
im = ax.imshow(abvalsarray, cmap=custom_cmap)

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(cases)), labels=cases)
ax.set_yticks(np.arange(len(gonames_list)), labels=gonames_list)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor", fontsize=8)

cbar = plt.colorbar(im)
cbar.ax.set_ylabel("<- 100% A compartment           100% B compartment ->",
                   rotation=270, horizontalalignment="center", labelpad=20)

# # Loop over data dimensions and create text annotations.
# for i in range(len(gonames_list)):
#     for j in range(len(cases)):
#         text = ax.text(j, i, abvalsarray[i, j],
#                        ha="center", va="center", color="w")

#fig.tight_layout()
plt.savefig(plot3d)
plt.clf()


# plot Fig4D ab heatmaps
# get list of goterms for heatmap by taking intersections of vennlists
combs = combinations(vennlistlist_spectral, 2)
intxlist = []
for i, j in combs:
    intx = list(set(i) & set(j))
    intxlist.extend(intx)
intxset = list(set(intxlist))  # final list of goterms to be used

print(len(intxset))
if len(intxset) > 30:
    intxset = intxset[0:29]
print(len(intxset))

# create ndarray for storing ab vals
arrayshape = (len(intxset), len(cases))
abvalsarray = np.ndarray(arrayshape)
# populate array with abvals
for gotermi in range(len(intxset)):
    for casei in range(len(cases)):
        try:
            print(cases[casei] + " - " + intxset[gotermi])
            abvalsarray[gotermi, casei] = goabreplsdictspectral[cases[casei]][intxset[gotermi]]  # pulls ab val from dict
        except KeyError as ke:
            print(ke)
            abvalsarray[gotermi, casei] = np.nan  # if for some reason val not found, print error and fill with no val

print(abvalsarray.T)

# add GO names to term list for use as axis labels
# make goid:goname dictionary
gonames_list = []
print("creating goid:goname dictionary...")
with open(goobo_file, "r") as gof:
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

for go in intxset:  # just use first dict for iterating through so order is consistent
    try:
        goname = termdict[go]  # fetch name of goterm from obo
        if len(goname) > 37:
            goname = goname[:38] + "..."  # truncate goname if it is too long so it fits on plot
        gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later
    except KeyError as ke:
        gonames_list.append(go + " - " + "Unknown")

# create heatmap
fig, ax = plt.subplots()
im = ax.imshow(abvalsarray, cmap=custom_cmap)

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(cases)), labels=cases)
ax.set_yticks(np.arange(len(gonames_list)), labels=gonames_list)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor", fontsize=8)

cbar = plt.colorbar(im)
cbar.ax.set_ylabel("<- 100% A compartment           100% B compartment ->",
                   rotation=270, horizontalalignment="center", labelpad=20)
# # Loop over data dimensions and create text annotations.
# for i in range(len(gonames_list)):
#     for j in range(len(cases)):
#         text = ax.text(j, i, abvalsarray[i, j],
#                        ha="center", va="center", color="w")

#fig.tight_layout()
plt.savefig(plot4d)
plt.clf()