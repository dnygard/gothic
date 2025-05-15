from gothic import *

for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
	gname = case + "FixedBins500kb_annotated.graphml"
	outname = case + "FixedBins500kb_COO.tsv"
	make_graph_coo(gname, outfile=outname)
