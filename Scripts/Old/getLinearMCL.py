from gothic import *
import ast


# compare average linear distances of nodes in same MCL/Spectral cluster versus average linear distance globally
def do_linear_analysis_mcl(graphfile, resultsfile, outfile, plotfile, gencodefile):

    # create "gene":"chr:startpos-endpos" dictionary
    print("creating gene:position dictionary...")
    geneposdict = {}  # initialize dict
    gencode = pd.read_csv(gencodefile)  # load gencode table
    for index, row in gencode.iterrows():  # iterate through table
        chrom = row['seqid']
        startpos = int(row['start'])
        endpos = int(row['end'])
        name = row['attributes'].split(";")[1].split("=")[1].split(".")[0]  # extracts ensembl gene id from attributes column

        valstring = chrom + ":" + str(startpos) + "-" + str(endpos)
        geneposdict[name] = valstring  # add gene to dict
    del gencode  # delete from memory to free up resources

    # filter geneposdict so only contains genes that are annotated in the network
    # load graph
    if type(graphfile) == str:
        g = load_graph(graphfile)
    elif type(graphfile) == Graph:
        g = graphfile
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
    # build graph genes list
    graphgeneslist = []
    for v in g.vertices():
        graphgeneslist.extend(g.vp.genes[v])
    graphgeneslist = list(set(graphgeneslist))
    print(graphgeneslist[1])
    print(list(geneposdict.keys())[1])
    print("# of genes in graph: " + str(len(graphgeneslist)))  # TODO remove later
    del g

    # use graphgeneslist to filter geneposdict
    toremovelist = []
    for gene in geneposdict.keys():
        if gene not in list(graphgeneslist):
            toremovelist.append(gene)
    for gene in toremovelist:
        geneposdict.pop(gene)

    # get lists of intersecting genes from gprofiler results file
    clustgenelist = []
    with open(resultsfile, "r") as rf:
        for line in rf:
            if not line.startswith(","):  # ignore header lines
                genesliststring = line.split("query")[1].split("]")[1].split("[")[1]  # get "intersections" line from gprofiler results
                geneslistraw = genesliststring.split(", ")  # hack off [] and split into list of strings
                geneslist = [x[1:-1] for x in geneslistraw]  # hack off '' to clean up gene names
                clustgenelist.append(geneslist)  # append list to list

    # for each cluster, get average linear distance for genes on same chromosome
    avg_distlist = []  # list of average distances for each enriched GOterm gene in each cluster
    for genelist in clustgenelist:
        distlist = []  # list of distances for averaging later
        # get all pairwise combinations
        vcombs = list(combinations(genelist, 2))

        # get distances for all gene pairs that exist on the same chromosome
        for g1, g2 in vcombs:
            if g1 in geneposdict and g2 in geneposdict:
                pos1 = geneposdict[g1]
                chr1 = pos1.split(":")[0]
                start1 = int(pos1.split(":")[1].split("-")[0])
                end1 = int(pos1.split(":")[1].split("-")[1])

                pos2 = geneposdict[g2]
                chr2 = pos2.split(":")[0]
                start2 = int(pos2.split(":")[1].split("-")[0])
                end2 = int(pos2.split(":")[1].split("-")[1])

                #print(pos1 + " -> " + pos2)  # TODO read as verbose option

                if chr1 == chr2:
                    if start1 > end2 and end1 > end2:  # if 1 is later in chrom than 2
                        dist = start1 - end2
                    elif start1 < start2 and end1 < start2:  # if 2 is later in chrom than 1
                        dist = start2 - end1
                    else:
                        dist = 0  # if any overlap dist is 0

                    distlist.append(dist)  # append distances
                else:
                    continue  # skip if not same chromosome
            else:
                print(g1 + " or " + g2 + " not in dict. Skipping...")  # TODO remove?
                continue

        # add avg dist to list
        if distlist:
            avgdist = sum(distlist)/len(distlist)  # get average (mean) of distances from all nodes that share chromosomes
        else:
            avgdist = -1  # if distlist is empty, no nodes are on same chrom so set avgdist to -1
        avg_distlist.append(avgdist)

    # get average linear distance for random pairs of genes (n samples = 10000)
    print("sampling random gene pairs to get average linear distance (N=10000)...")
    samplecounter = 0
    distlist = []
    n = 10000

    with tqdm(total=n) as pbar:  # make pbar for sampling loop
        while samplecounter < n:
            genepair = random.sample(list(geneposdict.keys()), 2)
            if geneposdict[genepair[0]].split(":")[0] == geneposdict[genepair[1]].split(":")[0]:  # if random genes share the same chromosome
                g1, g2 = genepair  # unpack tuple
                #print(genepair)
                #print(geneposdict[genepair[0]])
                #print(geneposdict[genepair[1]])
                start1 = geneposdict[genepair[0]].split(":")[1].split("-")[0]  # get start pos for gene 1
                stop1 = geneposdict[genepair[0]].split(":")[1].split("-")[1]  # get stop pos for gene 1
                start2 = geneposdict[genepair[1]].split(":")[1].split("-")[0]  # get start pos for gene 2
                stop2 = geneposdict[genepair[1]].split(":")[1].split("-")[1]  # get stop pos for gene 2
                if int(start1) > int(stop2):  # if gene 1 is later pos that gene 2, dist = start1 - stop2
                    pairdist = int(start1) - int(stop2)
                    #print(str(start1) + "-" + str(stop2))
                elif int(start2) > int(stop1):  # else is opposite
                    pairdist = int(start2) - int(stop1)
                    #print(str(start2) + "-" + str(stop1))
                else:  # only remaining case should be if start and stop are equal (which should be exceedingly rare)
                    pairdist = 0

                distlist.append(pairdist)  # append new distance to pairdistlist
                samplecounter = samplecounter + 1  # increment counter
                pbar.update()

    # write avg distances to file
    with open(outfile, "w") as outf:
        for dist in avg_distlist:
            outf.write(str(dist) + "\n")

    # plotting
    # get mean and standard deviation for random gene pairs
    mean = sum(distlist) / len(distlist)
    variance = sum([((x - mean) ** 2) for x in distlist]) / len(distlist)
    stdev = variance ** 0.5

    # create a plot of linear distance vs random ranking with line for avg linear distance and std dev
    # create lists to plot
    print("plotting...")
    ranking = range(len(avg_distlist))  # data for x-axis

    # Create a scatter plot
    plt.figure().set_figwidth(15)
    plt.scatter(ranking, avg_distlist, s=0.5, label='Linear Distance vs significant gene cluster', color='black')

    # Set labels for the x and y axes
    plt.xlabel('number of enriched gene (random)')
    plt.ylabel('Linear distance (bp)')

    # add lines for avg linear distance and standard deviations
    meanlabel = "mean = " + str(mean)
    upperlabel = "stddev = " + str(stdev)
    plt.axhline(y=mean, color='red', linestyle='-', label=meanlabel)  # mean
    plt.axhline(y=mean + stdev, color='red', linestyle='--', label=upperlabel)  # stddev upper bound
    plt.axhline(y=mean - stdev, color='red', linestyle='--')  # stddev lower bound

    # Show a legend
    plt.legend()

    # Display the plot
    plt.savefig(plotfile)
    print("plot saved as " + plotfile)


if __name__ == "__main__":
    for case in ["1a", "1b", "1c"]:
        for method in ["MCL", "Spectral"]:
            print(case + " " + method)
            graphfile = "MEC" + case + "TADgraphPt0175cut100kwindow_oneminusw.gt"
            resultsfile = "MEC" + case + "TADgraphPt0175cut100kwindow_" + method + "enrichments.csv"
            outfile = "MEC" + case + method + "_LinearVs3Dresults.txt"
            plotfile = "MEC" + case + method + "_LinearVs3Dplot.png"
            gencodefile = "gencode_pcgenes.csv"
            do_linear_analysis_mcl(graphfile, resultsfile, outfile, plotfile, gencodefile)
