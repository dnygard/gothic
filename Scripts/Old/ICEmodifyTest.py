from gothic import *

infile = "MEC1cFixedBins500kb_wmodified.gt"
outfile = "MEC1cFixedBins500kb_normedScaledModified.gt"
regplotfile = "MEC1cFixedBins500kb_normedScaledModified_EdgeWeights.png"
logplotfile = "MEC1cFixedBins500kb_normedScaledModified_EdgeWeightsLogScale.png"

# do modification
simple_invert_weights(infile, wannotation="raw_weight", outfile=outfile)

# plot first on normal scale then log scale
plot_edge_weight_distribution(outfile, binsize=100, outfile=regplotfile, log=False)
plot_edge_weight_distribution(outfile, binsize=100, outfile=logplotfile, log=True)
