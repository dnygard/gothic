# script for creating GO annotated Hi-C network and retrieving preliminary TPD clustering results

from gothic import *

"""
full run of the pipeline using Micro-C data for clustering analysis of Gene Ontology terms, 
from sam files to final results dataframe plus preliminary figures.
 
Includes results for WTPD, MCL, and Spectral-OPTICS clustering 

define a function to do the whole run, such that it can be picked up at any point if something fails
requires the following arguments, which carry through all steps:
runname: str name to use as a file prefix for all intermediate files and plots 
sam1: str name of first Hi-C sam file to be used for graph creation (default=None)
sam2: str name of second Hi-C sam file to be used for graph creation (default=None)(default=None)
filedir: directory containing all required files for running (listed below) (default=".")
binsize: int size in base pairs that each node in the network will represent if using fixed-width binning (default=80000)
vmax: int maximum number of nodes to consider in caclulating WTPD of GO terms (default=10000)
ncores: int number of cores to use in parallel processing of certain steps (default=1)
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

def microc_full_run(runname, sam1=None, sam2=None, filedir=".", binsize=80000, vmax=10000, ncores=1,
                    step=0, endstep=12, binfile=None, saveadjlist=True, filterlevel=0):
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
    graphname = runname + "_noannotation.graphml"
    mymap = filedir + "/HUMAN_9606_idmapping_selected.tsv"
    gencode = filedir + "/gencode_pcgenes.csv"
    goobo = filedir + "/go-basic.obo"
    reactomemap = filedir + "/Ensembl2Reactome_All_Levels.txt"
    newgraphname = runname + "_incompleteannotated.graphml"
    if step <= 3 <= endstep:
        print("Annotating graph...")
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

        if filterlevel > 0:
            print("filtering edges...")
            graphname = newgraphname
            fpstring = "Top" + str(int(filterlevel * 100)) + "p"
            newgraphname = runname + "_" + fpstring + "filtered.gt"
            create_top_edges_graph(graphname, threshold=filterlevel, outfile=newgraphname)
        else:
            g = load_graph(newgraphname)
            g.save(runname + ".gt")

        # print("Modifying edge weights...")
        # graphname = newgraphname
        # newgraphname = runname + "_filteredModifiedEdges.gt"
        # modify_weights(graphname, outfile=newgraphname)
        # stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(NormalizeStep): ', stop - start)  # Timer Block

    # calculate All Pairs Shortest Paths (APSP) and create distance matrix
    if step <= 5 <= endstep:
        if filterlevel > 0:
            newgraphname = runname + "_filtered.gt"
        else:
            newgraphname = runname + ".gt"

        distgraphname = runname + "_wdistances.gt"
        annotate_apsp(newgraphname, distgraphname)  # calculate all pairs shortest paths as graph annotations
        matfilename = runname + "_distmatrix.tsv"
        distmatrix_from_annotations(distgraphname, outfile=matfilename)  # convert/format distances to distance matrix
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(DijkstraStep): ', stop - start)  # Timer Block

    # plot edge weight and shortest paths distributions, save as base file name for all future uses
    if step <= 6 <= endstep:
        graphname = runname + "_wdistances.gt"
        newgraphname = runname + ".gt"  # naming convention assumes basic graph name is the complete one TODO is that ok?
        matfilename = runname + "_distmatrix.tsv"
        weightplotname = runname + "_edgeWeightPlot.png"
        apspplotname = runname + "_APSPdistributionPlot.png"

        print("Saving graph as " + newgraphname + "...")
        g = load_graph(graphname)
        g.save(newgraphname)
        del g

        print("Plotting edge weight and shortest paths distributions...")
        plot_edge_weight_distribution(graphname, outfile=weightplotname)
        plot_shortest_path_pairs_distribution(matfilename, outfile=apspplotname)

    # generate graph report
    if step <= 7 <= endstep:
        graphname = runname + ".gt"
        print("Generating graph report...")
        generate_graph_report(graphname, outfile=runname + "graphReport.txt")

    # get sampling vector and perform MC sampling
    if step <= 8 <= endstep:
        graphname = runname + ".gt"
        print("Starting MC sampling...")
        matfilename = runname + "_distmatrix.tsv"
        tvlengths = get_vlengths(graphname)  # gets list of k for mc sampling
        vlengths = [i for i in tvlengths if i <= vmax]  # trims list of k to only include vals less than vmax
        distrnfile = runname + "_distributions.csv"
        montecarlo_sample_tpds(matfilename, vlengths, graphname, ncores=ncores, outfile=distrnfile)
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(MCsamplingStep): ', stop - start)  # Timer Block

    # create shuffled graph
    if step <= 9 <= endstep:
        print("Shuffling graph annotations...")
        graphname = runname + ".gt"
        newgraphname = runname + "_goshuffled.gt"
        swap_all_goterms(graphname, outfile=newgraphname)
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(ShuffleStep): ', stop - start)  # Timer Block

    # get real and shuffled TPDs and p-values for all go terms
    if step <= 10 <= endstep:
        print("Getting real and shuffled TPDs and p-values...")
        graphname = runname + ".gt"
        newgraphname = runname + "_goshuffled.gt"
        matfilename = runname + "_distmatrix.tsv"
        gd = make_go_dict(graphname)  # need to make godict before doing all_go_tpd()
        shufgd = make_go_dict(newgraphname)  # need to make godict before doing all_go_tpd()
        all_go_tpd(gd, matfilename, runname + "_realTPDs.csv")
        all_go_tpd(shufgd, matfilename, runname + "_shuffledTPDs.csv")
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(TPDstep): ', stop - start)  # Timer Block;l''''''

    # get final dataframe with "nnodes", "tpd", "pval", "shuftpd", "shufpval" for all go terms
    if step <= 11 <= endstep:
        print("Collecting final results...")
        graphname = runname + ".gt"
        distrnfile = runname + "_distributions.csv"
        resultscsvname = runname + "_resultsDF.csv"
        topgolistname = runname + "_topGOterms.csv"
        fdrthresholdsname = runname + "_FDRthresholds.csv"

        df = get_go_tpd_pvals(graphname, runname + "_realTPDs.csv", runname + "_shuffledTPDs.csv", distrnfile)
        df.to_csv(resultscsvname)

        # also get list of top go terms with FDR
        get_top_goterms(resultscsvname, outfile=topgolistname)
        get_real_v_null_significants(resultscsvname, outfile=fdrthresholdsname)

        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(EndStep): ', stop - start)  # Timer Block

if __name__ == "__main__":

    gothic_full_run("YourRunName", sam1="YourSAM1.sam", sam2="YourSAM2.sam")