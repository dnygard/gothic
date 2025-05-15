from gothic import *


# takes a graph file and annotates it with
def label_compartments_no_gc(graph, outfile="labelled.gt"):
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







gfile = "MEC1cTADgraphPt0175cut100kwindow.gt"
cfile = "MEC1c_SpectralClusters.csv"
outg = "testab.gt"
#label_compartments_no_gc(gfile, outfile=outg)
ab_cluster_report(gfile, cfile)
