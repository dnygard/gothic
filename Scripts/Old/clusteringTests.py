from gothic import *
from sklearn.cluster import OPTICS
import markov_clustering as mcl
from scipy.sparse.linalg import eigsh


# performs spectral clustering on graph
def do_spectral_clustering(graph, outfile, method="OPTICS"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    # get graph spectrum for clustering
    adjacency_matrix = adjacency(g, weight=g.ep.weight).toarray()  # Compute the weighted adjacency matrix
    degree_matrix = np.diag(adjacency_matrix.sum(axis=0))  # Compute the degree matrix for weighted graph
    laplacian_matrix = degree_matrix - adjacency_matrix  # Compute the Laplacian matrix
    eigenvalues, eigenvectors = eigsh(laplacian_matrix, k=3, which='SM')  # Eigenvalue Decomposition

    # do OPTICS Clustering
    clust = OPTICS(min_samples=3, xi=0.05).fit(eigenvectors)  # do clustering
    # cluster_labels = clust.labels_[clust.ordering_]
    cluster_labels = clust.labels_

    clusternodedict = {}
    node = 0
    for i in [int(x) for x in cluster_labels]:
        if i == -1:
            node = node + 1
            continue
        else:
            if i in clusternodedict.keys():
                clusternodedict[i].append(node)
            else:
                clusternodedict[i] = [node]
            node = node + 1

    with open(outfile, "w") as f:
        for k in clusternodedict.keys():
            f.write(",".join([str(x) for x in clusternodedict[k]]) + "\n")


# performs OPTICS clustering
def do_optics_clustering(graph):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    # do OPTICS clustering on whole matrix
    adjacency_matrix = adjacency(g, weight=g.ep.raw_weight).toarray()  # Compute the weighted adjacency matrix
    clust = OPTICS(min_samples=3, xi=0.05).fit(adjacency_matrix)  # do clustering
    cluster_labels = clust.labels_[clust.ordering_]
    print(cluster_labels)

    return cluster_labels


# performs Markov clustering on graph
def do_markov_clustering(graph, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    matrix = adjacency(g, weight=g.ep.raw_weight).toarray()

    result = mcl.run_mcl(matrix)  # run MCL with default parameters
    clusters = mcl.get_clusters(result)  # get clusters

    for line in clusters:
        writeline = str(line).split("(")[-1].split(")")[0] + "\n"
        f.write(str(writeline))


# tests clusters from
def test_go_enrichment(x):
    print()


# creates upset plot from multiple lists of enriched go terms (i.e. comparing enrichment across cases)
def plot_enrichment_upsets():
    print()


if __name__ == "__main__":

    # for each case
    for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
        # for each filter level
        for flevel in [0.1, 0.2, 0.5]:
            # set file names
            fpstring = "Top" + str(int(flevel * 100)) + "p"
            run = case + "FixedBins500kb" + fpstring + "Edges"
            print(run)
            graphfile = case + "FixedBins500kb_edgeModifiedFiltered" + fpstring + ".gt"
            mclresultsfile = run + "_MCLclusters.csv"
            mcldistrnplotfile = run + "_MCLclusterDistrn.png"
            spectralresultsfile = run + "_SpectralClusters.csv"
            spectraldistrnplotfile = run + "_SpectralClusterDistrn.png"

            # # do mcl
            # clusters = do_markov_clustering(graphfile)
            # clusterhist = []
            # with open(mclresultsfile, "w") as f:
            #     for line in clusters:
            #         clusterhist.append(len(line))
            #         writeline = str(line).split("(")[-1].split(")")[0] + "\n"
            #         f.write(str(writeline))
            #
            #     # plot distrn
            #     plt.figure(figsize=(10, 6))
            #     plt.hist(clusterhist, bins=max(clusterhist)/2, color='black', alpha=0.8, rwidth=0.95)
            #     plt.xlabel('Size of cluster')
            #     plt.ylabel('Count')
            #     plt.savefig(mcldistrnplotfile, bbox_inches='tight')
            #     plt.clf()
            #     plt.close()
            #
            # print("mcl done")

            # do spectral + optics
            clusters = do_spectral_clustering(graphfile, spectralresultsfile)
            #np.savetxt(spectralresultsfile, clusters, delimiter=",")

            # plot distrn
            plt.figure(figsize=(10, 6))
            plt.hist(list(clusters), bins='auto', color='black', alpha=0.7, rwidth=0.85)
            plt.xlabel('Size of cluster')
            plt.ylabel('Count')
            plt.savefig(spectraldistrnplotfile, bbox_inches='tight')
            plt.clf()
            plt.close()

            print("spectral done")
