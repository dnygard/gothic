from gothic import *
from sklearn.cluster import OPTICS
import markov_clustering as mcl
from scipy.sparse.linalg import eigsh

# performs Markov clustering on graph
def do_markov_clustering(graph, inflation=2):
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

    result = mcl.run_mcl(matrix, inflation=inflation)  # run MCL with default parameters
    clusters = mcl.get_clusters(result)  # get clusters

    return clusters


if __name__ == "__main__":

    try:
        os.makedirs("inflationTests", mode=0o777, exist_ok=False)
    except FileExistsError as fe:
        print("directory exists, proceeding")

    for infl in [i / 10 for i in range(12, 51)]:
        print(infl)
        graphfile = "MEC1cFixedBins500kb_edgeModifiedFilteredTop10p.gt"
        mclresultsfile = "inflationTests/MEC1c500kbTop10p_MCLclusters_inflation" + str(int(infl*10)) + ".csv"
        mcldistrnplotfile = "inflationTests/MEC1c500kbTop10p_MCLclusterDistrn_inflation" + str(int(infl*10)) + ".png"
        # do mcl
        clusters = do_markov_clustering(graphfile, inflation=infl)
        clusterhist = []
        with open(mclresultsfile, "w") as f:
            for line in clusters:
                clusterhist.append(len(line))
                writeline = str(line).split("(")[-1].split(")")[0] + "\n"
                f.write(str(writeline))

            # plot distrn
            plt.figure(figsize=(10, 6))
            plt.hist(clusterhist, bins=int(max(clusterhist)/2), color='black', alpha=0.8, rwidth=0.95)
            plt.xlabel('Size of cluster')
            plt.ylabel('Count')
            plt.savefig(mcldistrnplotfile, bbox_inches='tight')
            plt.clf()
            plt.close()
