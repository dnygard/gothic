from itertools import combinations
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
from scipy.stats import pearsonr

# performs linear vs 3D clustering analysis
# first creates dictionary of "gene":"chr:startpos-endpos" dictionary
# plots linear distance of GO term genes vs TPD GO term ranking in terms of TPD
# plot also include line of true estimated average distance between any two genes

# for each go term
# get all nodes annotated with that term
# get all genes associated with that term from the nodes
# count number of pairs of genes that exist on the same chromosome (for each chromosome present)
# for each pair of genes calculate the linear distance between them
# find average distance for set by adding up all pairwise linear distances and dividing by the pair count


gencodefile = "gencode_pcgenes.csv"
resultsfile = "MEC2ndTADgraphPt0175cut100kwindow_topGOtermsWnames.csv"
gafname = "goa_human.gaf"
outname = "linearVs3DplotMEC2.png"

# create "gene":"chr:startpos-endpos" dictionary
print("creating gene:position dictionary...")
geneposdict = {}  # initialize dict
gencode = pd.read_csv(gencodefile)  # load gencode table
for index, row in gencode.iterrows():  # iterate through table
    chrom = row['seqid']
    startpos = int(row['start'])
    endpos = int(row['end'])
    name = row['attributes'].split(";")[3].split("=")[1]  # extracts gene name from attributes column

    valstring = chrom + ":" + str(startpos) + "-" + str(endpos)
    geneposdict[name] = valstring  # add gene to dict
del gencode  # delete from memory to free up resources

# create "term":[genes] dictionary from gaf file from GO consortium
print("creating GOterm:genes dictionary...")
gogenesdict = {}  # dict for storing "GOid":[list of genes] dictionary
with open(gafname, "r") as gaf:
    for line in gaf:
        if not line.startswith("!"):  # ignore header lines
            gene = line.split("\t")[2]  # retrieve gene name from 3rd column
            goid = line.split("\t")[4]  # retrieve GOid from 5th column
            if goid in gogenesdict:  # if go id is in dict already, append new gene to list
                gogenesdict[goid].append(gene)
            else:  # else create new dict entry for genes as list of strings (with only one entry for now)
                gogenesdict[goid] = [gene]

# get ranked list of go terms from results file
print(" getting ranked list of GO terms from results file...")
rankedgolist = []
with open(resultsfile, "r") as resf:  # open results file
    for line in resf:  # iterate through file to extract ranked list of go terms
        if not line.split(",")[0] == "goterm":  # skip header line of file
            rankedgolist.append(line.split(",")[0])  # append go term to list

# in order, iterate through ranked list of GO terms
print("getting linear distances for all pairs of genes per ranked GO terms...")
golindistdict = {}
for term in tqdm(rankedgolist):  # iterate through ranked list of goterms
    # get all genes associated with that term
    try:
        genelist = gogenesdict[term]
    except KeyError as ke:
        print("Could not get gene info for " + term)
        golindistdict[term] = -1
        continue
    # get position/loc info for each term in list
    locdict = {}  # dictionary will contain chromosomes as key, list of positions within that chromosome as value
    for g in genelist:
        try:
            locinfo = geneposdict[g]
        except KeyError as ke:
            print("KeyError: " + g + " not found in gencode dictionary")
            continue

        if locinfo.split(":")[0] in locdict:  # if chrom is already in dict, append new position info to list
            locdict[locinfo.split(":")[0]].append(locinfo.split(":")[1])
        else:  # else create new entry for chromosome
            locdict[locinfo.split(":")[0]] = [locinfo.split(":")[1]]

    # get number of valid gene pairs (pairs of genes that exist on the same chromosome)
    paircount = 0  # keep track of pair count so we can use as divisor for average pairwise distance
    dist = 0  # keep total distance for all pairs as single value, to be divided by pair count
    for chrom in locdict.keys():
        if len(locdict[chrom]) >= 2:  # only can count pairwise distance if there are at least two genes on chromosome
            uniquelist = list(set(locdict[chrom]))  # ensure all gene entries in list are unique by removing duplicates
            combs = list(combinations(uniquelist, 2))  # get all pairwise combinations
            paircount = paircount + len(combs)  # add pair count for chromosome to total pair count
            for c in combs:  # for each pair, find distance between them
                g1, g2 = c  # unpack tuple
                start1 = g1.split("-")[0]  # get start pos for gene 1
                stop1 = g1.split("-")[1]  # get stop pos for gene 1
                start2 = g2.split("-")[0]  # get start pos for gene 2
                stop2 = g2.split("-")[1]  # get stop pos for gene 2
                print("start1: " + str(int(start1)) + ", stop1: " + str(int(stop1)) + ", start2: " + str(int(start2))
                      + ", stop2: " + str(int(stop2)))

                if int(start1) > int(stop2):  # if gene 1 is later pos that gene 2, dist = start1 - stop2
                    pairdist = int(start1) - int(stop2)
                    print(str(start1) + "-" + str(stop2))
                elif int(start2) > int(stop1):  # else is opposite
                    pairdist = int(start2) - int(stop1)
                    print(str(start2) + "-" + str(stop1))
                else:  # only remaining case should be if start and stop are equal (which should be exceedingly rare)
                    pairdist = 0
                print("pairdist: " + str(pairdist))
                dist = dist + pairdist  # add distance of pair to total summed distance
            print(term + ": " + str(dist) + "bp")

    # get average linear distance for gene set
    try:
        avglindist = float(dist) / float(paircount)
        print("average linear distance for set: " + str(avglindist))
    except ZeroDivisionError as ze:
        print("Could not calculate linear distance for term: " + term)
        continue
    # add to dictionary of "GO":avglindist
    golindistdict[term] = avglindist

# get average linear distance for random pairs of genes (n samples = 1000000)
print("sampling random gene pairs to get average linear distance (N=1000000)...")
samplecounter = 0
distlist = []
n = 1000000
with tqdm(total=n) as pbar:  # make pbar for sampling loop
    while samplecounter < n:
        genepair = random.sample(list(geneposdict.keys()), 2)
        if geneposdict[genepair[0]].split(":")[0] == geneposdict[genepair[1]].split(":")[0]:  # if random genes share the same chromosome
            g1, g2 = genepair  # unpack tuple
            print(genepair)
            print(geneposdict[genepair[0]])
            print(geneposdict[genepair[1]])
            start1 = geneposdict[genepair[0]].split(":")[1].split("-")[0]  # get start pos for gene 1
            stop1 = geneposdict[genepair[0]].split(":")[1].split("-")[1]  # get stop pos for gene 1
            start2 = geneposdict[genepair[1]].split(":")[1].split("-")[0]  # get start pos for gene 2
            stop2 = geneposdict[genepair[1]].split(":")[1].split("-")[1]  # get stop pos for gene 2
            if int(start1) > int(stop2):  # if gene 1 is later pos that gene 2, dist = start1 - stop2
                pairdist = int(start1) - int(stop2)
                print(str(start1) + "-" + str(stop2))
            elif int(start2) > int(stop1):  # else is opposite
                pairdist = int(start2) - int(stop1)
                print(str(start2) + "-" + str(stop1))
            else:  # only remaining case should be if start and stop are equal (which should be exceedingly rare)
                pairdist = 0

            distlist.append(pairdist)  # append new distance to pairdistlist
            samplecounter = samplecounter + 1  # increment counter
            pbar.update()

# get mean and standard deviation for random gene pairs
mean = sum(distlist) / len(distlist)
variance = sum([((x - mean) ** 2) for x in distlist]) / len(distlist)
stdev = variance ** 0.5


# create a plot of linear distance vs GO term ranking with line for avg linear distance and std dev
# create lists to plot
print("plotting...")
ranking = range(len(rankedgolist))  # data for x-axis
rankeddistlist = []  # initialize list for y-axis
for term in rankedgolist:  # append distance for each ranked GO term in order from most to least significant
    try:
        rankeddistlist.append(golindistdict[term])
    except KeyError as ke:
        print("no linear distance found for " + term)
        rankeddistlist.append(-1)

# Create a scatter plot
plt.figure().set_figwidth(15)
plt.scatter(ranking, rankeddistlist, s=0.5, label='Linear Distance vs GO Ranking', color='black')

# Set labels for the x and y axes
plt.xlabel('Ranking of GO term by TPD significance')
plt.ylabel('Linear distance (bp)')

# add lines for avg linear distance and standard deviations
meanlabel = "mean = " + str(mean)
upperlabel = "stddev = " + str(stdev)
plt.axhline(y=mean, color='red', linestyle='-', label=meanlabel)  # mean
plt.axhline(y=mean+stdev, color='red', linestyle='--', label=upperlabel)  # stddev upper bound
plt.axhline(y=mean-stdev, color='red', linestyle='--')  # stddev lower bound

# Show a legend
plt.legend()

# Display the plot
plt.savefig(outname)
print("plot saved as " + outname)

corr, _ = pearsonr(ranking, rankeddistlist)

print("Pearson Correlation = " + str(corr))
