from gothic import *
from matplotlib.colors import LinearSegmentedColormap

cases = ["1a", "1b", "1c"]
plot2e = "MECTADgraphPt0175cut100kwindow_fig2eOLD.png"
plot3d = "MECTADgraphPt0175cut100kwindow_fig3d.png"
plot4d = "MECTADgraphPt0175cut100kwindow_fig4d.png"

#constants
fdrthresh = 0.2
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

    # #### Fig 2: TPD for all 3 Mammary Epithelial cell replicates
    #     # Figure 2: A) B) C) D) E)
    #
    # # read in results file as lists
    # with open(tpd_results_file, "r") as rf:
    #     gotermlist = []
    #     gonamelist = []
    #     nodeswithtermlist = []
    #     chromswithtermlist = []
    #     tpdlist = []
    #     pvallist = []
    #     fdrlist = []
    #     abratiolist = []
    #     avglineardistancelist = []
    #     shuftpdlist = []
    #     shufpvallist = []
    #
    #     for line in rf:
    #         if not line.startswith("GOterm"):  # skip header line
    #             splitline = line.strip().split(",")
    #             gotermlist.append(splitline[0])
    #             gonamelist.append(splitline[1])
    #             nodeswithtermlist.append(splitline[2])
    #             chromswithtermlist.append(splitline[3])
    #             tpdlist.append(float(splitline[4]))
    #             pvallist.append(float(splitline[5]))
    #             fdrlist.append(float(splitline[6]))
    #             abratiolist.append(float(splitline[7]))
    #             avglineardistancelist.append(float(splitline[8]))
    #             shuftpdlist.append(float(splitline[9]))
    #             shufpvallist.append(float(splitline[10]))
    #
    # # E) A/B compartments analysis with the ratio in a heatmap where all 3 replicates are combined (next to each other as columns) on the heatmap [ ]
    # # Store output in abtpdlistlist to be plotted after looping through all replicates
    # # need to store both GO names and associated ab ratios (tuple?)
    # for i in range(len(abratiolist)):
    #     goabreplsdicttpd[case][gotermlist[i]] = abratiolist[i]  # store ab values in dict of dicts
    #
    # # F) Venn Diagram of significant GO terms (need to choose an FDR based on above results) where each replicate is a circle. [ ]
    # # Supplemental 1: Revigo figures of significant GO terms (need to choose an FDR based on above results) [ ]
    #     # DONE VIA GUI, NOT DONE PROGRAMATICALLY HERE
    # # build list of lists to be sent to venn diagram software
    # vennlist = []
    # with open(tpd_vennlist, "w") as f:
    #     for i in range(len(fdrlist)):
    #         if fdrlist[i] <= fdrthresh:
    #             vennlist.append(gotermlist[i])
    #             f.write(gotermlist[i] + "\n")
    # vennlistlist.append(vennlist)
    #
    # # clear lists for tpd results from memory
    # del gotermlist
    # del gonamelist
    # del nodeswithtermlist
    # del chromswithtermlist
    # del tpdlist
    # del pvallist
    # del fdrlist
    # del abratiolist
    # del avglineardistancelist


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

    # # clear mcl lists from memory
    # del gotermlist
    # del gonamelist
    # del nodeswithterminclusterlist
    # del numberofnodeswithtermingraphlist
    # del numberofnodeswithterminclusterlist
    # del chromswithtermlist
    # del chromswithterminclusterlist
    # del pvallist
    # del abratiolist
    # del avglineardistancelist


#     # #### Fig 4: Spectral clustering and GO analysis for all 3 Mammary Epithelial cell replicates
#     # read in results file as lists
#     with open(spectral_results_file, "r") as rf:
#         gotermlist = []
#         gonamelist = []
#         nodeswithterminclusterlist = []
#         numberofnodeswithtermingraphlist = []
#         numberofnodeswithterminclusterlist = []
#         chromswithtermlist = []
#         chromswithterminclusterlist = []
#         pvallist = []
#         abratiolist = []
#         avglineardistancelist = []
#
#         for line in rf:
#             if not line.startswith("GOterm"):  # skip header line
#                 splitline = line.strip().split(",")
#                 gotermlist.append(splitline[0])
#                 nodeswithterminclusterlist.append(splitline[1])
#                 numberofnodeswithtermingraphlist.append(splitline[2])
#                 numberofnodeswithterminclusterlist.append(splitline[3])
#                 chromswithtermlist.append(splitline[4])
#                 chromswithterminclusterlist.append(splitline[5])
#                 pvallist.append(float(splitline[6]))
#                 abratiolist.append(float(splitline[7]))
#                 avglineardistancelist.append(float(splitline[8]))
#                 gonamelist.append(splitline[9])
#
#     # 4d) A/B compartments analysis with the ratio in a heatmap where all 3 replicates are combined (next to each other as
#     # columns) on the heatmap for significant GO terms. Ratio is calculated only on the set of nodes in a cluster annotated with the significant GO terms. [ ]
#     # DONE LATER SO ALL 3 REPLICATES CAN BE COMBINED
#     # Store output in abtpdlistlist to be plotted after looping through all replicates
#     # need to store both GO names and associated ab ratios (tuple?)
#     for i in range(len(abratiolist)):
#         goabreplsdictspectral[case][gotermlist[i]] = abratiolist[i]  # store ab values in dict of dicts
#
#     # 4e) Venn Diagram of significant GO terms (need to choose an FDR-adjusted enrichment p-value based on above results)
#     # where each replicate is a circle. [ ]
#
#     # build list of lists to be sent to venn diagram software
#     vennlist = []
#     with open(spectral_vennlist, "w") as f:
#         for i in range(len(pvallist)):
#             if pvallist[i] <= pvalthresh:
#                 vennlist.append(gotermlist[i])
#                 f.write(gotermlist[i] + "\n")
#     vennlistlist_spectral.append(vennlist)
#
#
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
#
# # plot 2E ab heatmaps
# # get list of goterms for heatmap by taking intersections of vennlists
# combs = combinations(vennlistlist, 2)
# intxlist = []
# for i, j in combs:
#     intx = list(set(i) & set(j))
#     intxlist.extend(intx)
# intxset = list(set(intxlist))  # final list of goterms to be used
#
# # create ndarray for storing ab vals
# arrayshape = (len(intxset), len(cases))
# abvalsarray = np.ndarray(arrayshape)
# # populate array with abvals
# for gotermi in range(len(intxset)):
#     for casei in range(len(cases)):
#         try:
#             abvalsarray[gotermi, casei] = goabreplsdicttpd[cases[casei]][intxset[gotermi]]  # pulls ab val from dict
#         except KeyError as ke:
#             print(ke)
#             abvalsarray[gotermi, casei] = np.nan  # if for some reason val not found, print error and fill with no val
#
# # add GO names to term list for use as axis labels
# # make goid:goname dictionary
# gonames_list = []
# print("creating goid:goname dictionary...")
# with open(goobo_file, "r") as gof:
#     termdict = {}  # dictionary for mapping go ids to terms (key=ID, value=term)
#     ids = []  # list to append ids and alt ids to
#     for line in gof:
#         splitline = line.split(": ")
#         if splitline[0] == "name":
#             goname = splitline[1].strip()
#         if splitline[0] == "id":
#             ids.append(splitline[1].strip())
#         if splitline[0] == "altid":
#             ids.append(splitline[1].strip())
#         if splitline[0] == "def":
#             for goid in ids:
#                 termdict[goid] = goname
#             ids = []  # empty list of ids
# for go in intxset:  # just use first dict for iterating through so order is consistent
#     goname = termdict[go]  # fetch name of goterm from obo
#     if len(goname) > 37:
#         goname = goname[:38] + "..."  # truncate goname if it is too long so it fits on plot
#     gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later
#
# # create heatmap
# fig, ax = plt.subplots()
# im = ax.imshow(abvalsarray)
#
# # Show all ticks and label them with the respective list entries
# ax.set_xticks(np.arange(len(cases)), labels=cases)
# ax.set_yticks(np.arange(len(gonames_list)), labels=gonames_list)
#
# # Rotate the tick labels and set their alignment.
# plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#          rotation_mode="anchor")
#
# cbar = plt.colorbar(im)
# cbar.ax.set_ylabel("<- 100% A compartment           100% B compartment ->",
#                    rotation=270, horizontalalignment="center", labelpad=20)
# # # Loop over data dimensions and create text annotations.
# # for i in range(len(gonames_list)):
# #     for j in range(len(cases)):
# #         text = ax.text(j, i, abvalsarray[i, j],
# #                        ha="center", va="center", color="w")
#
# fig.tight_layout()
# plt.savefig(plot2e)
# plt.clf()


# plot Fig3D ab heatmaps
# # get list of goterms for heatmap by taking intersections of vennlists
# combs = combinations(vennlistlist_mcl, 2)
# intxlist = []
# for i, j in combs:
#     intx = list(set(i) & set(j))
#     intxlist.extend(intx)
# intxset = list(set(intxlist))  # final list of goterms to be used
# # ^use this block if only doing heatmap for terms intersecting at least 2 cases

# # get list of goterms for heatmap by taking 3-way intersection of vennlists
# intxset = list(set(vennlistlist_mcl[0]) & set(vennlistlist_mcl[1]) & set(vennlistlist_mcl[2]))

# # get list of goterms for heatmap by taking intersection of only 1a and 1b vennlists, since 1c had no clusters
intxset = list(set(vennlistlist_mcl[0]) & set(vennlistlist_mcl[1]))

###  <- adjust set size here 3d
# adjust set size for plotting
numprint = 40  # non-hard-coded number of rows to print per heatmap

print(len(intxset))
ogintxset = intxset  # save orriginal intxset for comparing lengths during plotting loop
i1 = 0
i2 = numprint - 1
namenum = 1
while i2 < len(ogintxset) - 1:
    if len(ogintxset) > numprint:
        intxset = ogintxset[i1:i2]
    else:
        i2 = len(ogintxset) - 1
    print(len(intxset))

    # create ndarray for storing ab vals
    arrayshape = (len(intxset), len(cases)-1)
    abvalsarray = np.ndarray(arrayshape)
    # populate array with abvals
    for gotermi in range(len(intxset)):
        for casei in range(len(cases)-1):
            try:
                #print(cases[casei] + " - " + intxset[gotermi])
                abvalsarray[gotermi, casei] = goabreplsdictmcl[cases[casei]][intxset[gotermi]]  # pulls ab val from dict
            except KeyError as ke:
                #print(ke)
                abvalsarray[gotermi, casei] = np.nan  # if for some reason val not found, print error and fill with no val

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
            if len(goname) > 38:
                goname = goname[:37] + "..."  # truncate goname if it is too long so it fits on plot
            gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later
        except KeyError as ke:
            gonames_list.append(go + " - " + "Unknown")

    # create heatmap
    #plt.figure().set_figheight(100)
    fig, ax = plt.subplots()
    im = ax.imshow(abvalsarray)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(cases)-1), labels=cases[0:2])
    ax.set_yticks(np.arange(len(gonames_list)), labels=gonames_list, fontsize=6)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor", fontsize=6)

    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel("<- 100% A compartment           100% B compartment ->",
                       rotation=270, horizontalalignment="center", labelpad=20)

    # # Loop over data dimensions and create text annotations.
    # for i in range(len(gonames_list)):
    #     for j in range(len(cases)):
    #         text = ax.text(j, i, abvalsarray[i, j],
    #                        ha="center", va="center", color="w")

    #fig.tight_layout()
    newplotname = plot3d.split(".")[0] + "_plot" + str(namenum) + "OLD." + plot3d.split(".")[1]
    plt.savefig(newplotname)
    plt.clf()

    i1 = i1 + numprint
    if i2 + numprint <= len(ogintxset) - 1:
        i2 = i2 + numprint
    else:
        i2 = len(ogintxset) - 1
    namenum = namenum + 1

# # plot Fig4D ab heatmaps
# # get list of goterms for heatmap by taking intersections of vennlists
# combs = combinations(vennlistlist_spectral, 2)
# intxlist = []
# for i, j in combs:
#     intx = list(set(i) & set(j))
#     intxlist.extend(intx)
# intxset = list(set(intxlist))  # final list of goterms to be used
# ^use this block if only doing heatmap for terms intersecting at least 2 cases

# # get list of goterms for heatmap by taking 3-way intersection of vennlists
# intxset = list(set(vennlistlist_spectral[0]) & set(vennlistlist_spectral[1]) & set(vennlistlist_spectral[2]))
#
# ###  <- adjust set size here 4d
# # adjust set size for plotting
# numprint = 40  # non-hard-coded number of rows to print per heatmap
#
# print(len(intxset))
# ogintxset = intxset  # save orriginal intxset for comparing lengths during plotting loop
# i1 = 0
# i2 = numprint - 1
# namenum = 1
# while i2 < len(ogintxset) - 1:
#     if len(ogintxset) > numprint:
#         intxset = ogintxset[i1:i2]
#     else:
#         i2 = len(ogintxset) - 1
#     print(len(intxset))
#
#     # create ndarray for storing ab vals
#     arrayshape = (len(intxset), len(cases))
#     abvalsarray = np.ndarray(arrayshape)
#     # populate array with abvals
#     for gotermi in range(len(intxset)):
#         for casei in range(len(cases)):
#             try:
#                 #print(cases[casei] + " - " + intxset[gotermi])
#                 abvalsarray[gotermi, casei] = goabreplsdictspectral[cases[casei]][intxset[gotermi]]  # pulls ab val from dict
#             except KeyError as ke:
#                 #print(ke)
#                 abvalsarray[gotermi, casei] = np.nan  # if for some reason val not found, print error and fill with no val
#
#     # add GO names to term list for use as axis labels
#     # make goid:goname dictionary
#     gonames_list = []
#     print("creating goid:goname dictionary...")
#     with open(goobo_file, "r") as gof:
#         termdict = {}  # dictionary for mapping go ids to terms (key=ID, value=term)
#         ids = []  # list to append ids and alt ids to
#         for line in gof:
#             splitline = line.split(": ")
#             if splitline[0] == "name":
#                 goname = splitline[1].strip()
#             if splitline[0] == "id":
#                 ids.append(splitline[1].strip())
#             if splitline[0] == "altid":
#                 ids.append(splitline[1].strip())
#             if splitline[0] == "def":
#                 for goid in ids:
#                     termdict[goid] = goname
#                 ids = []  # empty list of ids
#
#     for go in intxset:  # just use first dict for iterating through so order is consistent
#         try:
#             goname = termdict[go]  # fetch name of goterm from obo
#             if len(goname) > 38:
#                 goname = goname[:37] + "..."  # truncate goname if it is too long so it fits on plot
#             gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later
#         except KeyError as ke:
#             gonames_list.append(go + " - " + "Unknown")
#
#     # create heatmap
#     #plt.figure().set_figheight(100)
#     fig, ax = plt.subplots()
#     im = ax.imshow(abvalsarray)
#
#     # Show all ticks and label them with the respective list entries
#     ax.set_xticks(np.arange(len(cases)), labels=cases)
#     ax.set_yticks(np.arange(len(gonames_list)), labels=gonames_list, fontsize=6)
#
#     # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#              rotation_mode="anchor", fontsize=6)
#
#     cbar = plt.colorbar(im)
#     cbar.ax.set_ylabel("<- 100% A compartment           100% B compartment ->",
#                        rotation=270, horizontalalignment="center", labelpad=20)
#
#     # # Loop over data dimensions and create text annotations.
#     # for i in range(len(gonames_list)):
#     #     for j in range(len(cases)):
#     #         text = ax.text(j, i, abvalsarray[i, j],
#     #                        ha="center", va="center", color="w")
#
#     #fig.tight_layout()
#     newplotname = plot4d.split(".")[0] + "_plot" + str(namenum) + "OLD." + plot4d.split(".")[1]
#     plt.savefig(newplotname)
#     plt.clf()
#
#     i1 = i1 + numprint
#     if i2 + numprint <= len(ogintxset) - 1:
#         i2 = i2 + numprint
#     else:
#         i2 = len(ogintxset) - 1
#     namenum = namenum + 1



