from gothic import *

# parse genome fasta according to bins in binbedfile, then count and output G/C ratios for each bin
# the regions file better be sorted in ascending order by chromosome and chrom position or else
def get_node_gc_ratios(binbedfile, genomefastafile, outfile=None):

    bf = open(binbedfile, "r")  #
    gf = open(genomefastafile, "r")  # me and who?

    # iterate through bedfile to get lists of node start positions
    bin_chrom_start_list = []

    for line in bf:
        bin_chrom_start_list.append(int(line.split()[1]))
    bf.close()

    # iterate through genome file and count GC content for each node, at each bin boundary store the GC vs ATGC ratio in a list
    node_gc_ratios_list = []  # list for storing the gc content ratio (NUMGC/NUMALLBASES) of each region
    pos_counter = 0
    base_counter = 0
    gc_counter = 0
    done_flag = 0
    nextpos = bin_chrom_start_list.pop(0)

    for line in gf:

        if done_flag:  # if no more bins to get gc ratios for, break loop
            break

        if line.startswith(">"):
            if base_counter != 0:
                gc_ratio = float(gc_counter)/float(base_counter)
                node_gc_ratios_list.append(gc_ratio)
                pos_counter = 0
                base_counter = 0
                gc_counter = 0

        else:
            base_list = [x for x in line]
            base_list.pop()  # remove newline character from line list so it is not counted as a base
            for b in base_list:
                # add base to counts depending on identity
                base_counter = base_counter + 1
                if b == "G" or b == "C":
                    gc_counter = gc_counter + 1

                pos_counter = pos_counter + 1

                if pos_counter >= nextpos:
                    gc_ratio = float(gc_counter) / float(base_counter)
                    node_gc_ratios_list.append(gc_ratio)
                    #print(str(gc_counter) + "/" + str(base_counter))  #TODO REMOVE WHEN DONE
                    #print(gc_ratio)  #TODO REMOVE WHEN DONE

                    base_counter = 0
                    gc_counter = 0

                    if bin_chrom_start_list:  # get next startpos if there is one
                        nextpos = bin_chrom_start_list.pop(0)
                        #print(nextpos)  #TODO REMOVE WHEN DONE

                    else:  # else raise done flag and break loop
                        done_flag = 1
                        break

    gf.close()

    if outfile:
        with open(outfile, "w") as outf:
            for gc in node_gc_ratios_list:
                outf.write(str(gc) + "\n")

    return node_gc_ratios_list


# annotate graph nodes with G/C ratio as a vertex property "gc_ratio"
def annotate_nodes_with_gc(graph, gclist, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    gc_vp = g.new_vertex_property("double")

    for i in range(0, len(gclist)):
        gc_vp[g.vertex(i)] = gclist[i]

    g.vp["gc_ratio"] = gc_vp
    g.save(outfile)


# takes a graph file and annotates it with value of first principal component "ab_eigenval" for getting A/B compartments
def annotate_compartments_no_gc(graph, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    adjacency_matrix = adjacency(g, weight=g.ep.weight).toarray()  # Compute the weighted adjacency matrix
    pca = PCA(n_components=1)
    X = pca.fit_transform(adjacency_matrix)

    print(X)
    print(type(X))
    print(X.shape)
    xlist = [float(x[0]) for x in X.tolist()]
    print(xlist)
    print(len(xlist))

    ab_vp = g.new_vertex_property("double")

    for i in range(len(xlist)):
        ab_vp[g.vertex(i)] = xlist[i]

    g.vp["ab_eigenval"] = ab_vp
    g.save(outfile)


# only use if annotated with above function
def ab_cluster_report(graph, clusterfile, outfile="clusters_ab_report.tsv"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    with open(outfile, "w") as outf:
        outf.write("numV\tnumNeg\tnumPos\tposRatio\n")
        with open(clusterfile, "r") as f:
            for line in f:  # for each node in cluster
                vlist = line.split(",")
                ablist = []

                for v in vlist:  # get corresponding eigenvalues of these nodes
                    ablist.append(g.vp.ab_eigenval[v])

                numv = str(len(ablist))
                numpos = str(len([x for x in ablist if x >= 0]))
                numneg = str(len([x for x in ablist if x < 0]))
                posratio = str(float(numpos)/float(numv))

                outf.write(numv + "\t" + numneg + "\t" + numpos + "\t" + posratio + "\n")


# uses ab_eigenval and gc_ratio vertex properties to assign A/B compartment labels to nodes as vp "ab_compartment"
def determine_ab_compartments(graph, outfile, verbose=False):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    if verbose:
        # get number of nodes for use in progress bars
        vcount = 0
        for v in g.vertices():
            vcount = vcount+1

    # build list of nodes for neg and pos compartments
    neg_vlist = []
    pos_vlist = []
    if verbose:
        print("Building compartment nodelists...")
        pbar = tqdm(total=vcount)
    for v in g.vertices():
        if g.vp.ab_eigenval[v] >= 0:
            pos_vlist.append(v)
        else:
            neg_vlist.append(v)
        if verbose:
            pbar.update()

    # get gc ratios for each of those lists
    neg_gclist = []
    pos_gclist = []

    # determine which compartment is A and which is B
    if verbose:
        print("Determining compartment gc ratios...")
    for v in neg_vlist:  # get list of GC ratios for negative compartment
        neg_gclist.append(g.vp.gc_ratio[v])
    for v in pos_vlist:  # get list of GC ratios for positive compartment
        pos_gclist.append(g.vp.gc_ratio[v])

    # get the average GC content ratio for each compartment
    neg_gc_avg = sum(neg_gclist) / len(neg_gclist)
    if verbose:
        print("Negative compartment avg GC ratio is " + str(neg_gc_avg) + "...")
    pos_gc_avg = sum(pos_gclist) / len(pos_gclist)
    if verbose:
        print("Positive compartment avg GC ratio is " + str(pos_gc_avg) + "...")

    # assign "A" to whichever commpartment has the higher GC ratio, "B" to other
    if neg_gc_avg > pos_gc_avg:
        neg = "A"  # negative values are A compartment
        pos = "B"  # positive values are B compartment
        if verbose:
            print("Negative is A, Positive is B...")

    elif pos_gc_avg > neg_gc_avg:
        neg = "B"  # negative values are B compartment
        pos = "A"  # positive values are A compartment
        if verbose:
            print("Positive is A, Negative is B...")
    else:
        print("Could not determine A/B compartments. Stopping...")
        sys.exit()

    # create new vertex property "ab_compartment" and fill it with compartment labels
    ab_vp = g.new_vertex_property("string")

    if verbose:
        print("Assigning A/B compartment status to nodes...")
        pbar = tqdm(total=vcount)
    for v in g.vertices():
        if g.vp.ab_eigenval[v] >= 0:  # label positive eigenvalues with compartment determined above
            ab_vp[v] = pos
        else:  # label negative eigenvalues with compartment determined above
            ab_vp[v] = neg
        if verbose:
            pbar.update()

    # make internal
    g.vp.ab_compartment = ab_vp

    # save graph
    if verbose:
        print("Saving...")
    g.save(outfile)


if __name__ == "__main__":

    genome_file = "GRCh38_noalt_as.fa"
    cases = ["1a", "1b", "1c"]
    go_ab_repls_dict = {}  # dict of {"goname":{"1a":ABval}, {"1b":ABval}, {"1c":ABval}}
    go_ab_array_list = []  # list of lists for later conversion to array via np.column_stack()
    gonames_list = []  # list of gonames for later use as heatmap labels
    go_list_file = "MEC1a1b1cSharedGOsTPD.txt"
    obofile = "go-basic.obo"

    # make heatmap
    # read in list of go terms to plot (overlap of multiple replicates in this case)
    go_list = []
    with open(go_list_file, "r") as f:
        for line in f:
            go_list.append(line.strip())
    print(go_list)  # TODO remove

    # make goid:goname dictionary
    print("creating goid:goname dictionary...")
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
    ## check that worked
    #for key in termdict.keys():
    #    print(key + "->" + termdict[key])
    #print(len(termdict.keys()))
    #print(termdict.keys())


    for case in cases:
        graph_file = "MEC" + case + "TADgraphPt0175cut100kwindow.gt"
        bin_file = "TADsMEC" + case + "0.0175cut100000window.tads"
        top_terms_file = "MEC" + case + "TADgraphPt0175cut100kwindow_topGOtermsPLUSNAMESPLUSAB.csv"
        print(case)  # TODO remove

        # # annotate graph with A/B compartments
        # print("Getting A/B compartments for " + case + "...")
        # print("Getting GC content ratios...")
        # gc_list = get_node_gc_ratios(bin_file, genome_file)
        # print("Annotating graph with pseudocompartments and gc ratios...")
        # annotate_nodes_with_gc(graph_file, gc_list, graph_file)
        # annotate_compartments_no_gc(graph_file, graph_file)
        # print("Annotating with true A/B compartments...")
        # determine_ab_compartments(graph_file, graph_file, verbose=True)
        # print("Done " + case + "!")

        # get vals for go_ab_repls_dict
        print("fetching AB values...")
        go_ab_repls_dict[case] = {}  # dict of {"goname":{"1a":ABval}, {"1b":ABval}, {"1c":ABval}}
        g = load_graph(graph_file)
        for go in go_list:
            vlist = get_nodes_with_term(g, go)
            ablist = [g.vp.ab_compartment[i] for i in vlist]
            abratio = sum(i == "A" for i in ablist) / len(ablist)
            go_ab_repls_dict[case][go] = abratio

        #print(go_ab_repls_dict[case])

    # convert to array for heatmap and list of "GO:GONAME" strings to use as names
    print("evoking AB array...")
    for go in go_ab_repls_dict[list(go_ab_repls_dict.keys())[0]].keys():  # just use first dict for iterating through so order is consistent
        abs_list = []
        goname = termdict[go]  # fetch name of goterm from obo
        gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later
        for repl in go_ab_repls_dict.keys():  # for each case get the AB value at that go, and add it to a list
            abs_list.append(go_ab_repls_dict[repl][go])

        # add list of vals for this go term to list of lists to later
        go_ab_array_list.append(abs_list)

    # convert go_ab_array_list to array
    print("plotting...")
    heatmap_vals_array = np.column_stack(go_ab_array_list)
    np.savetxt("heatmap_vals_copy.csv", heatmap_vals_array, delimiter=",")  # save for posterity

    # create heatmap
    vegetables = gonames_list

    farmers = cases

    harvest = heatmap_vals_array

    fig, ax = plt.subplots()
    im = ax.imshow(harvest)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(farmers)), labels=farmers)
    ax.set_yticks(np.arange(len(vegetables)), labels=vegetables)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(vegetables)):
        for j in range(len(farmers)):
            text = ax.text(j, i, harvest[i, j],
                           ha="center", va="center", color="w")

    fig.tight_layout()
    plt.savefig("MEC_AB_heatmap_test.png")
    plt.show()




    # create heatmap where row labels are a list of clustered go terms, columns are replicates, colour is A/B compartment ratio

