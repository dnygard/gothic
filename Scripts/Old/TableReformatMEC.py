from gothic import *

# generate/reformat tables from TPD, MCL, Spectral experiments for figure generation and publication

if __name__ == '__main__':
    cases = '1a', '1b', '1c'
    for case in cases:
        print(case)  # TODO remove?
        # instantiate all files f
        graphfile = "MEC" + case + "TADgraphPt0175cut100kwindow_oneminusw.gt"
        tpd_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_topGOterms.csv"
        tpd_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_SupplementalTable.csv"
        tpd_og_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_resultsDF.csv"
        tpd_realvshuf_file = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_SignificantClustersRealVsShuffled.csv"
        tpd_update_file = "MEC" + case + "TADgraphPt0175cut100kwindow_TPD_SupplementalTableFIXED.csv"
        mcl_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLenrichments.csv"
        mcl_cluster_file = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters.csv"
        mcl_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCLclusters_Min3Filtered.csv"
        mcl_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_SupplementalTable.csv"
        mcl_plotfile = "MEC" + case + "TADgraphPt0175cut100kwindow_MCL_LinDistPlot.png"
        spectral_results_file = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectralenrichments.csv"
        spectral_cluster_file = "MEC" + case + "TADgraphPt0175cut100kwindow_SpectralClusters.csv"
        spectral_cluster_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_SpectralClusters_Min3Filtered.csv"
        spectral_outfile = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_SupplementalTable.csv"
        spectral_plotfile = "MEC" + case + "TADgraphPt0175cut100kwindow_Spectral_LinDistPlot.png"
        gaf_file = "goa_human.gaf"
        gencode_file = "gencode_pcgenes.csv"
        golist_file = "MEC" + case + "TADgraphPt0175cut100kwindow_AnnotatedGOterms.txt"
        goobo_file = "go-basic.obo"

        # ## filter mcl/spectral cluster files (REMOVED IN FAVOUR OF FILTERING PER CLUSTER AT RUNTIME)
        # # define filtering function
        # def filter_mcl_clusterfile(clusterfile, outfile, minsize=3):
        #     with open(clusterfile, "r") as cf:
        #         with open(outfile, "w") as outf:
        #             for line in cf:
        #                 if len(line.split(",")) < minsize:
        #                     continue
        #                 else:
        #                     outf.write(line)
        # # execute filtering
        # filter_mcl_clusterfile(mcl_cluster_file, mcl_cluster_outfile)
        # mcl_cluster_file = mcl_cluster_outfile
        minclustersize = 3

        ## load graph
        g = load_graph(graphfile)

        ## get list of all GO terms in the network
        print("getting annotated GO terms...")  # TODO remove?
        golist = []
        for v in g.vertices():
            # unpack GO vertex annotation and add all terms to list
            vgolist = g.vp.goterms[v]
            for gt in vgolist:
                if gt not in golist:
                    golist.append(gt)
        # save golist to file
        with open(golist_file, "w") as f:
            for go in golist:
                f.write(go + "\n")

        ## parse go obo file to get goterm:name dict
        gonamedict = {}
        print("parsing GO obo file...")  # TODO remove?
        with open(goobo_file, "r") as gf:
            goid = "NULL"
            goname = "NULL"
            for line in gf:
                if line.startswith("[T"):
                    gonamedict[goid] = goname
                elif line.startswith("id:"):
                    goid = line.split(": ")[1].strip()
                elif line.startswith("name:"):
                    goname = line.split(": ")[1].strip()

        # # check dict
        # for k in list(gonamedict.keys()):
        #     print(k + ": " + gonamedict[k])

        ## get linear distances to append afterwards
        tpd_lindists = get_linear_distances_tpd(tpd_results_file, gaf_file, gencode_file, nsamples=10000, verbose=False,
                                                plot=False)

        # get FDRs from tpd_realvshuf_file
        print("Getting shuftpd and shufpval frrom file...")
        with open(tpd_realvshuf_file, "r") as rf:
            fdrreflist = []  # list to reference for fdr at each pval
            print("retrieving pval -> fdr list from tpd_realvshuf_file")
            for line in tqdm(rf):
                if not line.startswith("p"):  # skip header line
                    fdrreflist.append(line.split(", ")[3])
        fdrreflist.reverse()  # reverse list so order is from 0 -> 1
        print(len(fdrreflist))

        ### TPD get lists for...
        print("Getting TPD results for " + case + "...")  # TODO remove?
        writelist = ["GOterm", "NodesWithTerm", "ChromsWithTerm", "TPD", "Pval", "FDR",
                      "ABratio", "AvgLinearDistance", "ShufTPD", "ShufPval", "GOname"]  # create list of list for later file writing
        writeline = ",".join(writelist) + "\n"
        # get shuftpd, shufpval dict from tpd_og_results_file
        # create {goterm:shufresults} dictionary
        df = pd.read_csv(tpd_og_results_file)
        df = df.transpose()
        df.columns = df.iloc[0]
        df["goterm"] = df.index
        df = df[1:]
        df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                        "shufpval": "float", "goterm": "str"})
        goshuftpddict = {}  # dictionary to store shufpval and shuftpd for each go term, which acts as the key
        goshufpvaldict = {}  # dictionary to store shufpval and shuftpd for each go term, which acts as the key
        for index, row in df.iterrows():
            goshufpvaldict[index] = str(row["shufpval"])
            goshuftpddict[index] = str(row["shuftpd"])

        with open(tpd_results_file, "r") as rf:
            with open(tpd_outfile, "w") as of:
                of.write(writeline)
                fdrindex = 0
                for line in tqdm(rf):
                    if not line.startswith(","):  # skip header line
                        splitline = line.split(",")
                        # GO term
                        goterm = splitline[1]
                        # GO name
                        try:
                            goname = gonamedict[goterm].replace(",", "")
                        except KeyError as ke:
                            print("GO term not found in dict")
                            print(ke)
                            goname = "Unknown"

                        # Number of nodes with GO annotation in network
                        gocount = 0
                        chromlist = []
                        gonodelist = []
                        for v in g.vertices():
                            vgolist = g.vp.goterms[v]  # get list of goterms for vertex
                            if goterm in vgolist:
                                gocount = gocount + 1  # gocount is what we will use for results
                                chromlist.append(g.vp.vname[v].split(":")[0])  # get chromosome from vname
                                gonodelist.append(int(v))  # add vertex with annotation to list

                        # participating nodes
                        gonodelist = [str(x) for x in gonodelist]
                        strgonodelist = " ".join(gonodelist)

                        # Number of chromosomes with GO annotation in network
                        chroms = " ".join(list(set(chromlist)))

                        # TPD
                        TPD = line.split(",")[2]

                        # p-value
                        pval = line.split(",")[3]

                        # use fdrindex to get linear distance
                        thislindist = tpd_lindists[fdrindex]

                        # FDR
                        # thisfdr = newfdrlist[fdrindex]  # use this val for fdr
                        # fdrindex = fdrindex + 1  # use this val for fdr
                        ## get fdr from fdrreflist
                        fdrstep = 0.0001
                        num_digits = 4  # num digits to round to (should match fdrstep)
                        fdrindex = int(float(pval) * (1/fdrstep))
                        print("index: " + str(fdrindex))
                        thisfdr = fdrreflist[fdrindex]

                        # A/B ratio
                        # get A/B compartment ratio of all nodes in cluster
                        vlist = get_nodes_with_term(graphfile, goterm)
                        g = load_graph(graphfile)
                        ablist = [g.vp.ab_compartment[i] for i in vlist]
                        abratio = sum(i == "A" for i in ablist) / len(ablist)

                        # shuftpd
                        shuftpd = goshuftpddict[goterm]

                        # shufpval
                        shufpval = goshufpvaldict[goterm]

                        # append all stats to writelist
                        writelist = [goterm, strgonodelist, chroms, str(TPD), str(pval), str(thisfdr),
                                    str(abratio), str(thislindist), shuftpd, shufpval, goname]
                        writeline = ",".join(writelist) + "\n"
                        print(writeline)  # TODO remove
                        of.write(writeline)  # write writelist to file
                        of.flush()

            ## do monotonic transformation for FDRs
            with open(tpd_results_file, "r") as rf:
                fdrlist = []
                for line in rf:
                    if not line.startswith(","):
                        fdrlist.append(line.split(",")[4])
            newfdrlist = []
            bestval = 1
            fdrlist.reverse()  # so values are descending
            fdrlist = [float(x) for x in fdrlist]  # coerce str vals in list to float
            for fdr in fdrlist:
                if fdr < bestval:
                    bestval = fdr
                    newfdrlist.append(bestval)
                else:
                    newfdrlist.append(bestval)
            newfdrlist.reverse()  # now newfdrlist contains monotonically tfed vals in proper order


        ### MCL get lists for...

        print("Getting MCL results...")  # TODO remove?
        # parse cluster file as list of lists
        clustlist = [[0, 0, 0]]  # initialize clustlist with dummy nodelist at index 0
        with open(mcl_cluster_file, "r") as cf:
            for line in cf:
                clustlist.append(line.strip().split(","))  # at each index is a list of nodes belonging to that cluster

        # linear distance of the nodes with the GO term in the cluster
        mcl_lindists = do_linear_analysis_mcl(graphfile, mcl_results_file, gencode_file, nsamples=10000, verbose=True,
                                              plot=True, plotfile=mcl_plotfile)

        writelist = ["GOterm", "NodesInCluster", "NumberOfNodesWithTermInGraph",
                      "NumberOfNodesWithTermInCluster", "ChromsWithTerm", "ChromsWithTermInCluster", "Pval", "ABratio",
                      "AvgLinDistOfGenesInCluster", "GOname"]
        writeline = ",".join(writelist) + "\n"

        with open(mcl_outfile, "w+") as f:
            f.write(writeline)  # write header line
            clustcount = 0
            with open(mcl_results_file, "r") as rf:
                for line in rf:
                    if line.startswith(","):  # new header line means is a new cluster
                        clustcount = clustcount + 1
                        continue  # skip all other steps and move to next line in results file
                    # also skip line if current cluster size is less than minsize
                    if len(clustlist[clustcount]) < minclustersize:
                        continue

                    # parse line from results file and extract important columns

                    # GO term
                    goterm = line.split(",")[2]

                    # GO name (could just retrieve from results file, but gonames use commas which is annoying)
                    try:
                        goname = gonamedict[goterm].replace(",", "")
                    except KeyError as ke:
                        print("GO term not found in dict")
                        print(ke)
                        goname = "Unknown"

                    # Cluster ID = clustcount

                    # Number of nodes in the cluster
                    nclustnodes = len(clustlist[clustcount])

                    # Number of nodes with GO annotation in the network
                    gocount = 0
                    chromlist = []
                    gonodelist = []
                    for v in g.vertices():
                        vgolist = g.vp.goterms[v]  # get list of goterms for vertex
                        if goterm in vgolist:
                            gocount = gocount + 1  # gocount is what we will use for results
                            chromlist.append(g.vp.vname[v].split(":")[0])  # get chromosome from vname
                            gonodelist.append(int(v))  # add vertex with annotation to list

                    ngonodes = len(gonodelist)

                    # Number of nodes with GO annotation in the cluster
                    gocount = 0
                    clustchromlist = []
                    gonodelist = []
                    for v in clustlist[clustcount]:
                        vgolist = g.vp.goterms[v]  # get list of goterms for vertex
                        if goterm in vgolist:
                            gocount = gocount + 1  # gocount is what we will use for results
                            clustchromlist.append(g.vp.vname[v].split(":")[0])  # get chromosome from vname
                            gonodelist.append(int(v))  # add vertex with annotation to list

                    ngoclustnodes = len(gonodelist)

                    # Number of chromosomes with GO annotation in network
                    chroms = " ".join(list(set(chromlist)))

                    # Number of chromosomes with GO annotation in cluster
                    clustchroms = " ".join(list(set(clustchromlist)))

                    # Enrichment p-value
                    gotp = False
                    ibuffer = 0
                    while not gotp:

                        try:  # check if value at this index is a p-value (numeric)
                            float(line.split(",")[4 + ibuffer])  # check current index
                            pval = line.split(",")[4 + ibuffer]  # assign pval
                            gotp = True  # raise gotp flag
                        except ValueError as ve:
                            print(ve)
                            ibuffer = ibuffer + 1  # if not numeric, try after the next comma (next index)
                            if ibuffer > 100:  # do buffer max as break case
                                pval = np.nan
                                break
                            else:
                                continue


                    # FDR-adjusted enrichment p-value
                    # ^ This is the only pval given. As far as I can tell there is no way to see uncorrected gprofiler pvals

                    # A/B ratio of the nodes with the GO term in the cluster
                    # get A/B compartment ratio of all nodes in cluster
                    # vlist = get_nodes_with_term(graphfile, goterm)
                    # clustgointxlist = list(set(vlist) & set(clustlist[clustcount]))  # get intersection of clust and goterm vlists
                    # ablist = [g.vp.ab_compartment[i] for i in clustgointxlist]
                    # if len(ablist) > 0:
                    #     abratio = sum(i == "A" for i in ablist) / len(ablist)
                    # else:
                    #     abratio = -1

                    # ^ instead of this, just do ab ratio of whole cluster
                    ablist = [g.vp.ab_compartment[i] for i in clustlist[clustcount]]
                    abratio = sum(i == "A" for i in ablist) / len(ablist)

                    # add all vals to writelist
                    writelist = [goterm, nclustnodes, ngonodes, ngoclustnodes, chroms, clustchroms,
                                      pval, abratio, mcl_lindists[clustcount], goname]
                    writelist = [str(x) for x in writelist]
                    writeline = ",".join(writelist) + "\n"
                    f.write(writeline)  # write line to file
                    f.flush()


        ### Spectral get lists for...
        print("Getting spectral results...")  # TODO remove?
        # parse cluster file as list of lists
        clustlist = [[0, 0, 0]]  # initialize clustlist with dummy nodelist at index 0
        with open(spectral_cluster_file, "r") as cf:
            for line in cf:
                clustlist.append(
                    line.strip().split(","))  # at each index is a list of nodes belonging to that cluster

        # linear distance of the nodes with the GO term in the cluster
        spectral_lindists = do_linear_analysis_mcl(graphfile, spectral_results_file, gencode_file, nsamples=10000,
                                                   verbose=True, plot=True, plotfile=spectral_plotfile)

        writelist = ["GOterm", "NodesInCluster", "NumberOfNodesWithTermInGraph",
                     "NumberOfNodesWithTermInCluster", "ChromsWithTerm", "ChromsWithTermInCluster", "Pval", "ABratio",
                     "AvgLinDistOfGenesInCluster", "GOname"]
        writeline = ",".join(writelist) + "\n"
        with open(spectral_outfile, "w+") as f:
            f.write(writeline)
            clustcount = 0
            with open(spectral_results_file, "r") as rf:
                for line in rf:
                    if line.startswith(","):  # new header line means is a new cluster
                        clustcount = clustcount + 1
                        continue  # skip all other steps and move to next line in results file
                    # also skip line if current cluster size is less than minsize
                    if len(clustlist[clustcount]) < minclustersize:
                        continue

                    # parse line from results file and extract important columns

                    # GO term
                    goterm = line.split(",")[2]

                    # GO name
                    try:
                        goname = gonamedict[goterm].replace(",", "")
                    except KeyError as ke:
                        print("GO term not found in dict")
                        print(ke)
                        goname = "Unknown"

                    # Cluster ID = clustcount

                    # Number of nodes in the cluster
                    nclustnodes = len(clustlist[clustcount])

                    # Number of nodes with GO annotation in the network
                    gocount = 0
                    chromlist = []
                    gonodelist = []
                    for v in g.vertices():
                        vgolist = g.vp.goterms[v]  # get list of goterms for vertex
                        if goterm in vgolist:
                            gocount = gocount + 1  # gocount is what we will use for results
                            chromlist.append(g.vp.vname[v].split(":")[0])  # get chromosome from vname
                            gonodelist.append(int(v))  # add vertex with annotation to list

                    ngonodes = len(gonodelist)

                    # Number of nodes with GO annotation in the cluster
                    gocount = 0
                    clustchromlist = []
                    gonodelist = []
                    for v in clustlist[clustcount]:
                        vgolist = g.vp.goterms[v]  # get list of goterms for vertex
                        if goterm in vgolist:
                            gocount = gocount + 1  # gocount is what we will use for results
                            clustchromlist.append(g.vp.vname[v].split(":")[0])  # get chromosome from vname
                            gonodelist.append(int(v))  # add vertex with annotation to list

                    ngoclustnodes = len(gonodelist)

                    # Number of chromosomes with GO annotation in network
                    chroms = " ".join(list(set(chromlist)))

                    # Number of chromosomes with GO annotation in cluster
                    clustchroms = " ".join(list(set(clustchromlist)))

                    # Enrichment p-value
                    gotp = False
                    ibuffer = 0
                    while not gotp:

                        try:  # check if value at this index is a p-value (numeric)
                            float(line.split(",")[4 + ibuffer])  # check current index
                            pval = line.split(",")[4 + ibuffer]  # assign pval
                            gotp = True  # raise gotp flag
                        except ValueError:
                            ibuffer = ibuffer + 1  # if not numeric, try after the next comma (next index)
                            if ibuffer > 100:  # do buffer max as break case
                                print(ve)
                                pval = np.nan
                                break
                            else:
                                continue


                    # A/B ratio of the nodes with the GO term in the cluster
                    # get A/B compartment ratio of all nodes in cluster
                    # vlist = get_nodes_with_term(graphfile, goterm)
                    # clustgointxlist = list(set(vlist) & set(clustlist[clustcount]))  # get intersection of clust and goterm vlists
                    # ablist = [g.vp.ab_compartment[i] for i in clustgointxlist]
                    # if len(ablist) > 0:
                    #     abratio = sum(i == "A" for i in ablist) / len(ablist)
                    # else:
                    #     abratio = -1

                    # ^ instead of this, just do ab ratio of whole cluster
                    ablist = [g.vp.ab_compartment[i] for i in clustlist[clustcount]]
                    abratio = sum(i == "A" for i in ablist) / len(ablist)

                    # add all vals to writelist
                    writelist = [goterm, nclustnodes, ngonodes, ngoclustnodes, chroms, clustchroms,
                                      pval, abratio, spectral_lindists[clustcount], goname]
                    writelist = [str(x) for x in writelist]
                    writeline = ",".join(writelist) + "\n"
                    f.write(writeline)  # write line to file
                    f.flush()
