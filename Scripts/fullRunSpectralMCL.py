# script for creating GO annotated Hi-C network and retrieving preliminary MCL and Spectral-OPTICS clustering results

from gothic import *

"""
full run of the pipeline usig both Markov clustering and Spectral clustering with OPTICS for clustering analysis of 
Gene Ontology terms, from sam files to final results dataframe plus preliminary figures

define a function to do the whole run, such that it can be picked up at any point if something fails
requires the following arguments, which carry through all steps:
runname: str name to use as a file prefix for all intermediate files and plots 
sam1: str name of first Hi-C sam file to be used for graph creation (default=None)
sam2: str name of second Hi-C sam file to be used for graph creation (default=None)(default=None)
filedir: directory containing all required files for running (listed below) (default=".")
binsize: int size in base pairs that each node in the network will represent if using fixed-width binning (default=80000)
step: int step in pipeline to start the run at (default=0)
endstep: int step in pipeline to end the run at (default=12)
binfile: str name of regions file to use for node binning if using TAD bins or custom regions (default=None)
saveadjlist: bool if True, save the adjacency list as an intermediate file in case run fails (default=True)
filterlevel: int percent of top weighted edges to keep after network creation (default=0)

# files that need to be in filedir:
# - map(HUMAN_9606_idmapping_selected.tsv)
# - gencode(gencode_pcgenes.csv)
# - go obo(go-basic.obo)
# - Ensembl2Reactome_All_Levels.txt (mapping file only required if doing REACTOME annotation rather than GO)

TO RUN: enter your desired arguments into the run script at the end of this file, then execute as a python script
"""

def gothic_full_clustering_run(runname, sam1=None, sam2=None, filedir=".", binsize=80000,
                                step=0, endstep=12, binfile=None, saveadjlist=False, filterlevel=0.5):
    start = timeit.default_timer()

    if binfile:
        bindict = bins_to_dict(binfile)
        print("Bin dict created from file...")
    else:
        bindict = None

    # create graph
    if step <= 1 <= endstep:
        print("Creating Graph...")
        adjlist = sam_adjacencies(sam1, sam2, m=binsize, bins=bindict, verbose=True)
        print("Adjacency list created...")
        if saveadjlist:  # write adjacency list to file if saveadjlist option is true
            print("Writing adjacency list to file...")
            adjlistoutfile = runname + "_ADJLIST.tsv"
            save_adjlist(adjlist, adjlistoutfile)

    # convert adjlist to graphml
    if step <= 2 <= endstep:
        print("Converting adjlist to .graphml...")
        graphname = runname + "_noannotation.graphml"
        if saveadjlist:
            adjlist = []  # lil hack to overwrite the adjlist or else create one if there is none (ie start w this step)
            del adjlist  # clear up memory from filled adjlist dictionary object
            adjlistoutfile = runname + "_ADJLIST.tsv"
            adjlist_file_to_graphml(adjlistoutfile, outfile=graphname)
        else:
            hic_adjlist_to_graphml(adjlist, fileprefix=graphname)
            del adjlist  # clear up memory from filled adjlist dictionary object
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(GraphStep): ', stop - start)  # Timer Block

    # annotate graph
    if step <= 3 <= endstep:
        print("Annotating graph...")
        graphname = runname + "_noannotation.graphml"
        mymap = filedir + "/HUMAN_9606_idmapping_selected.tsv"
        gencode = filedir + "/gencode_pcgenes.csv"
        goobo = filedir + "/go-basic.obo"
        reactomemap = filedir + "/Ensembl2Reactome_All_Levels.txt"
        newgraphname = runname + "_incompleteannotated.graphml"
        genes_go_annotate(graphname, mapfile=mymap, gencodefile=gencode, m=binsize,
                          binfile=binfile, outfile=newgraphname, go_obo=goobo)
        oldgraphname = newgraphname
        newgraphname = runname + "_annotated.graphml"
        reactome_annotate(oldgraphname, mapfile=reactomemap, outfile=newgraphname)

        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(AnnotateStep): ', stop - start)  # Timer Block

    # normalize, adjust weights, filter, extract largest component
    if step <= 4 <= endstep:
        print("Normalizing...")
        graphname = runname + "_incompleteannotated.graphml"
        newgraphname = runname + "_normalized.graphml"
        cooname = runname + "prenormCOO.txt"
        normcooname = runname + "normedCOO.txt"
        prenorm_filter(graphname, graphname)  #filter edges with only one contact before normalizing
        make_graph_coo(graphname, outfile=cooname)
        ice_balance(cooname, normcooname)
        update_weights_from_coo(graphname, normcooname, newgraphname)

        print("filtering edges...")
        graphname = newgraphname
        fpstring = "Top" + str(int(filterlevel * 100)) + "p"
        newgraphname = runname + "_" + fpstring + "filtered.gt"
        create_top_edges_graph(graphname, threshold=filterlevel, outfile=newgraphname)
        print('Global Runtime(NormalizeStep): ', stop - start)  # Timer Block

    # plot edge weights and make graph report
    if step <= 5 <= endstep:
        graphname = newgraphname
        weightplotname = runname + "_edgeWeightPlot.png"
        generate_graph_report(graphname, outfile=runname + "graphReport.txt")
        plot_edge_weight_distribution(graphname, outfile=weightplotname)

    # do spectral optics and mcl clustering
    if step <= 6 <= endstep:
        graphname = newgraphname
        fpstring = "Top" + str(int(filterlevel * 100)) + "p"
        spectralresultsfile = runname + "_" + fpstring + "SpectralClusters.csv"
        mclresultsfile = runname + "_MCLclusters.csv"

        # do clustering
        do_markov_clustering(graphname, mclresultsfile)
        do_spectral_clustering(graphname, spectralresultsfile)

        # do clustering enrichment
        mcl_plotfile = runname + "_" + fpstring + "MCLenrichments.csv"
        spectral_plotfile = runname + "_" + fpstring + "Spectralenrichments.csv"

        get_go_enrichments(graphname, mclresultsfile, mcl_plotfile)
        get_go_enrichments(graphname, spectralresultsfile, spectral_plotfile)

if __name__ == "__main__":

    gothic_full_clustering_run("YourRunName", sam1="YourSAM1.sam", sam2="YourSAM2.sam")