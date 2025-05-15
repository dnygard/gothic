from gothic import *

# performs Markov clustering on graph
def do_markov_test(graph):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    matrix = adjacency(g, weight=g.ep.weight).toarray()

    # perform clustering using different inflation values from 1.5 and 2.5
    # for each clustering run, calculate the modularity
    for inflation in [i / 10 for i in range(11, 15)]:
        result = mcl.run_mcl(matrix, inflation=inflation)  # run MCL
        clusters = mcl.get_clusters(result)  # get clusters
        #Q = mcl.modularity(matrix=result, clusters=clusters)
        print(inflation)
        print(clusters)
        #print("inflation:", inflation, "modularity:", Q)

do_markov_test("MEC1aTADgraphPt0175cut100kwindow_oneminusw.gt")