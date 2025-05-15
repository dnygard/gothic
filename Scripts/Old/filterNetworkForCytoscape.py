from gothic import *

graphfile = "MEC1aTADgraphPt0175cut100kwindow_oneminusw.graphml"
outgraph = "MEC1aTADgraphPt0175cut100kwindow_oneminusw_top025percentGraph.graphml"
reportfile = "MEC1aTADgraphPt0175cut100kwindow_oneminusw_top025percentGraph_graphReport.txt"

create_top_edges_graph(graphfile, threshold=0.00025, outfile=outgraph, weights_are_distances=False)
generate_graph_report(outgraph, outfile=reportfile)
