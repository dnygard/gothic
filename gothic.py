import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import csv
import sys
from graph_tool.all import *
import concurrent.futures
import timeit
from tqdm import tqdm
import random
import math
import plotly.graph_objects as go
import time
from goatools import obo_parser
import subprocess
from scipy import stats as sps
import datetime
import os
import importlib
import seaborn as sns
import fanc
from scipy.sparse import coo_matrix, coo_array, csr_matrix, csc_matrix, vstack
from scipy.stats import gaussian_kde, ks_2samp, pearsonr
from scipy.sparse.linalg import eigsh
from iced import normalization
from gprofiler import GProfiler
from sklearn.cluster import OPTICS
from sklearn.decomposition import PCA
import markov_clustering as mcl
import hicstraw
import math
import logging

if sys.version_info[0] < 3:
    import StringIO
else:
    from io import StringIO

# GOTHiC: Gene Ontology Topology from Hi-C data

# utility function to quickly get line count of file for tracking progress through file iteration
def count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


# Converts Hi-C SAM files to adjacency list with nodes representing bins of desired size
# to use assemblies other than Hg38.13, define a custom vector of size 24 where entries are chrom lengths in order below
# to use custom bins, define a dict with keys chr1, chr2, ..., and values are a list of bin end positions
def sam_adjacencies(sam1, sam2, m=40000, bins=None, chromlengths=None, verbose=False):  # samt stands for "SAM Table", from the output of samTable function

    # Chromosome lengths based on Hg38.13
    if chromlengths is None:
        cr1 = 248956422
        cr2 = 242193529
        cr3 = 198295559
        cr4 = 190214555
        cr5 = 181538259
        cr6 = 170805979
        cr7 = 159345973
        cr8 = 145138636
        cr9 = 138394717
        cr10 = 133797422
        cr11 = 135086622
        cr12 = 133275309
        cr13 = 114364328
        cr14 = 107043718
        cr15 = 101991189
        cr16 = 90338345
        cr17 = 83257441
        cr18 = 80373285
        cr19 = 58617616
        cr20 = 64444167
        cr21 = 46709983
        cr22 = 50818468
        crX = 156040895
        crY = 57227415

    else:
        cr1 = chromlengths[0]
        cr2 = chromlengths[1]
        cr3 = chromlengths[2]
        cr4 = chromlengths[3]
        cr5 = chromlengths[4]
        cr6 = chromlengths[5]
        cr7 = chromlengths[6]
        cr8 = chromlengths[7]
        cr9 = chromlengths[8]
        cr10 = chromlengths[9]
        cr11 = chromlengths[10]
        cr12 = chromlengths[11]
        cr13 = chromlengths[12]
        cr14 = chromlengths[13]
        cr15 = chromlengths[14]
        cr16 = chromlengths[15]
        cr17 = chromlengths[16]
        cr18 = chromlengths[17]
        cr19 = chromlengths[18]
        cr20 = chromlengths[19]
        cr21 = chromlengths[20]
        cr22 = chromlengths[21]
        crX = chromlengths[22]
        crY = chromlengths[23]

    # chromosome list for checking whether RNAME value is valid or not
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    # initialize empty lists of lists of lists
    # list1: chromosomes (length 24: 22 + 23(X) +24(Y)), list2: bins(length = length(chr)/m),list3:weighted adjacencies
    # list sizes assume Hg38.13
    if bins is None:
        adjlist = [None] * 24
        adjlist[0] = [None] * int((cr1 / m) + 1)
        adjlist[1] = [None] * int((cr2 / m) + 1)
        adjlist[2] = [None] * int((cr3 / m) + 1)
        adjlist[3] = [None] * int((cr4 / m) + 1)
        adjlist[4] = [None] * int((cr5 / m) + 1)
        adjlist[5] = [None] * int((cr6 / m) + 1)
        adjlist[6] = [None] * int((cr7 / m) + 1)
        adjlist[7] = [None] * int((cr8 / m) + 1)
        adjlist[8] = [None] * int((cr9 / m) + 1)
        adjlist[9] = [None] * int((cr10 / m) + 1)
        adjlist[10] = [None] * int((cr11 / m) + 1)
        adjlist[11] = [None] * int((cr12 / m) + 1)
        adjlist[12] = [None] * int((cr13 / m) + 1)
        adjlist[13] = [None] * int((cr14 / m) + 1)
        adjlist[14] = [None] * int((cr15 / m) + 1)
        adjlist[15] = [None] * int((cr16 / m) + 1)
        adjlist[16] = [None] * int((cr17 / m) + 1)
        adjlist[17] = [None] * int((cr18 / m) + 1)
        adjlist[18] = [None] * int((cr19 / m) + 1)
        adjlist[19] = [None] * int((cr20 / m) + 1)
        adjlist[20] = [None] * int((cr21 / m) + 1)
        adjlist[21] = [None] * int((cr22 / m) + 1)
        adjlist[22] = [None] * int((crX / m) + 1)
        adjlist[23] = [None] * int((crY / m) + 1)

    else:
        adjlist = [None] * 24
        adjlist[0] = [None] * len(bins['chr1'])
        adjlist[1] = [None] * len(bins['chr2'])
        adjlist[2] = [None] * len(bins['chr3'])
        adjlist[3] = [None] * len(bins['chr4'])
        adjlist[4] = [None] * len(bins['chr5'])
        adjlist[5] = [None] * len(bins['chr6'])
        adjlist[6] = [None] * len(bins['chr7'])
        adjlist[7] = [None] * len(bins['chr8'])
        adjlist[8] = [None] * len(bins['chr9'])
        adjlist[9] = [None] * len(bins['chr10'])
        adjlist[10] = [None] * len(bins['chr11'])
        adjlist[11] = [None] * len(bins['chr12'])
        adjlist[12] = [None] * len(bins['chr13'])
        adjlist[13] = [None] * len(bins['chr14'])
        adjlist[14] = [None] * len(bins['chr15'])
        adjlist[15] = [None] * len(bins['chr16'])
        adjlist[16] = [None] * len(bins['chr17'])
        adjlist[17] = [None] * len(bins['chr18'])
        adjlist[18] = [None] * len(bins['chr19'])
        adjlist[19] = [None] * len(bins['chr20'])
        adjlist[20] = [None] * len(bins['chr21'])
        adjlist[21] = [None] * len(bins['chr22'])
        adjlist[22] = [None] * len(bins['chrX'])
        adjlist[23] = [None] * len(bins['chrY'])

    if verbose:
        # quickly get file line count for progress bar
        with open(sam1, 'rb') as fp:
            c_generator = count_generator(fp.raw.read)
            # count each \n
            count = sum(buffer.count(b'\n') for buffer in c_generator)
            print('SAM file line count:', count + 1)

    # iterate through SAM rows, and extract RNAME (chromosome) and POS, which is converted to a bin assignment
    with open(sam1) as file1, open(sam2) as file2:

        if verbose:
            pbar = tqdm(total=count+1)  # if verbose option set to True, create a progress bar that tracks parsed lines

        for line1, line2 in zip(file1, file2):
            try:
                if not line1.startswith('@'):  # ignores header lines of SAM file
                    # split input line by whitespace
                    s1 = line1.split()
                    s2 = line2.split()

                    # get info for first contact
                    chrom1 = s1[2]  # extracts chromosome name from 'RNAME' column of SAM
                    pos1 = int(s1[3])  # extracts leftmost mapping position from 'POS' column of SA
                    if chrom1 == '*':  # * means RNAME is unknown, so skip these lines
                        continue

                    bin1 = -1  # initialize bin1 as -1, then reassign according to binning logic
                    if bins is None:  # by default, bin is found by dividing the position by bin size
                        bin1 = int(pos1 / m)
                    else:  # otherwise loop through the bins to find which one it belongs in (by bin index)
                        bincounter = 0
                        for whichbin in bins[chrom1]:
                            binstart, binend = whichbin.split("-")
                            binstart = int(binstart)
                            binend = int(binend)
                            if binstart <= pos1 <= binend:
                                bin1 = bincounter
                                break  # exit loop when correct bin is found
                            bincounter = bincounter + 1
                        if bin1 == -1:  # if bin1 is not reassigned to a bin, exit this loop and go to next line
                            continue

                    # get info for second contact
                    chrom2 = s2[2]  # extracts chromosome name from 'RNAME' column of SAM
                    pos2 = int(s2[3])  # extracts leftmost mapping position from 'POS' column of SAM
                    if chrom2 == '*':  # * means RNAME is unknown, so skip these lines
                        continue

                    bin2 = -1
                    if bins is None:  # by default, bin is found by dividing the position by bin size
                        bin2 = int(pos2 / m)
                    else:  # otherwise loop through the bins to find which one it belongs in (by bin index)
                        bincounter = 0
                        for whichbin in bins[chrom2]:
                            binstart, binend = whichbin.split("-")
                            binstart = int(binstart)
                            binend = int(binend)
                            if binstart <= pos2 <= binend:
                                bin2 = bincounter
                                break  # exit loop when correct bin is found
                            bincounter = bincounter + 1
                        if bin2 == -1:  # if bin2 is not reassigned to a bin, exit this loop and go to next line
                            continue

                    if not (bin1 == bin2 and chrom1 == chrom2):  # filter out reads assigned to same bin
                        if any(item == chrom1 for item in chrlist) and any(item == chrom2 for item in chrlist):
                            # elif block that appends new connection/updates weight for appropriate bins
                            if chrom1 == 'chr1':
                                adjlist[0][bin1] = check_contact(adjlist[0][bin1], chrom2, bin2)
                            elif chrom1 == 'chr2':
                                adjlist[1][bin1] = check_contact(adjlist[1][bin1], chrom2, bin2)
                            elif chrom1 == 'chr3':
                                adjlist[2][bin1] = check_contact(adjlist[2][bin1], chrom2, bin2)
                            elif chrom1 == 'chr4':
                                adjlist[3][bin1] = check_contact(adjlist[3][bin1], chrom2, bin2)
                            elif chrom1 == 'chr5':
                                adjlist[4][bin1] = check_contact(adjlist[4][bin1], chrom2, bin2)
                            elif chrom1 == 'chr6':
                                adjlist[5][bin1] = check_contact(adjlist[5][bin1], chrom2, bin2)
                            elif chrom1 == 'chr7':
                                adjlist[6][bin1] = check_contact(adjlist[6][bin1], chrom2, bin2)
                            elif chrom1 == 'chr8':
                                adjlist[7][bin1] = check_contact(adjlist[7][bin1], chrom2, bin2)
                            elif chrom1 == 'chr9':
                                adjlist[8][bin1] = check_contact(adjlist[8][bin1], chrom2, bin2)
                            elif chrom1 == 'chr10':
                                adjlist[9][bin1] = check_contact(adjlist[9][bin1], chrom2, bin2)
                            elif chrom1 == 'chr11':
                                adjlist[10][bin1] = check_contact(adjlist[10][bin1], chrom2, bin2)
                            elif chrom1 == 'chr12':
                                adjlist[11][bin1] = check_contact(adjlist[11][bin1], chrom2, bin2)
                            elif chrom1 == 'chr13':
                                adjlist[12][bin1] = check_contact(adjlist[12][bin1], chrom2, bin2)
                            elif chrom1 == 'chr14':
                                adjlist[13][bin1] = check_contact(adjlist[13][bin1], chrom2, bin2)
                            elif chrom1 == 'chr15':
                                adjlist[14][bin1] = check_contact(adjlist[14][bin1], chrom2, bin2)
                            elif chrom1 == 'chr16':
                                adjlist[15][bin1] = check_contact(adjlist[15][bin1], chrom2, bin2)
                            elif chrom1 == 'chr17':
                                adjlist[16][bin1] = check_contact(adjlist[16][bin1], chrom2, bin2)
                            elif chrom1 == 'chr18':
                                adjlist[17][bin1] = check_contact(adjlist[17][bin1], chrom2, bin2)
                            elif chrom1 == 'chr19':
                                adjlist[18][bin1] = check_contact(adjlist[18][bin1], chrom2, bin2)
                            elif chrom1 == 'chr20':
                                adjlist[19][bin1] = check_contact(adjlist[19][bin1], chrom2, bin2)
                            elif chrom1 == 'chr21':
                                adjlist[20][bin1] = check_contact(adjlist[20][bin1], chrom2, bin2)
                            elif chrom1 == 'chr22':
                                adjlist[21][bin1] = check_contact(adjlist[21][bin1], chrom2, bin2)
                            elif chrom1 == 'chrX':
                                adjlist[22][bin1] = check_contact(adjlist[22][bin1], chrom2, bin2)
                            elif chrom1 == 'chrY':
                                adjlist[23][bin1] = check_contact(adjlist[23][bin1], chrom2, bin2)
                            else:
                                pass

            except KeyError:  # if mapping is for non-canonical chromosome, skip it
                continue

            if verbose:
                pbar.update(1)  # update pbar (if there is one) after each line

    if verbose:
        pbar.close()  # close pbar if one was created

    return adjlist


# converts adjacency list from sam_adjacencies() to graphml file
def hic_adjlist_to_graphml(adjlist, fileprefix="HicGraph"):
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    mlname = ""
    mlname = mlname.join([fileprefix, ".graphml"])

    g = Graph(directed=False)

    # create adjacency list to then import into graph with add_edgelist()
    contactlist = []  # for storing contacts
    nodelist = []  # for storing  unique node names, to later be  enumerated
    # initialize graph properties
    vname = g.new_vertex_property("string")  # vp map for node names (chromosome-bin pairs)
    weight = g.new_edge_property("double")  # ep map for edge weight

    for i in range(len(adjlist)):  # for each chromosome
        for j in range(len(adjlist[i])):  # for each bin
            # coerce chromosome-index pair to create node string (node u)
            nodeu = ""
            nodeu = nodeu.join([chrlist[i], ":", str(j)])  # get name of first node

            if adjlist[i][j]:
                for k in range(len(adjlist[i][j])):
                    # coerce chromosome-index pair to create node string (node v)
                    nodev = ""
                    nodev = nodev.join([adjlist[i][j][k][0], ":", str(adjlist[i][j][k][1])])  # get name of 2nd node

                    nodelist.append(nodeu)
                    nodelist.append(nodev)

                    # add edge between node string from outer loop and each node from here
                    contactlist.append([nodeu, nodev, weight])

    # once nodes and edges have been established, enumerate list of nodes and create graph using (index1,index2,weight)
    nodelist = list(set(nodelist))  # keep only unique entries in list of nodenames
    nlist = list(enumerate(nodelist))  # enumerate list for passing of index to graph-tool with add_edge_list()
    nodeindexmap = {}  # empty dict for mapping value:index pairs
    g.add_vertex(n=len(nodelist))  # add all vertices to graph first

    for n in nlist:
        nodeindexmap[n[1]] = n[0]  # add node index:name mapping
        vname[n[0]] = n[1]          # add name to vertex properties

    # prune edge list so a->b, b->a duplicates are combined and their weights are summed
    # uses a dictionary for ease of ooking up/ modifying weight value
    nameweightdict = {}
    for i in contactlist:
        n1, n2, w = i  # unpack "tuple" (actually a list rn but)
        # initialize strings for forward and reverse edge direction
        namestring = ""
        namestring = namestring.join([n1, "-", n2])
        rvrstring = ""
        rvrstring = namestring.join([n2, "-", n1])

        # if one name is already in the dictionary, update weight instead of adding
        if namestring in nameweightdict or rvrstring in nameweightdict:
            if namestring in nameweightdict:
                nameweightdict[namestring] = nameweightdict[namestring] + w
            if rvrstring in nameweightdict:
                nameweightdict[rvrstring] = nameweightdict[rvrstring] + w
        # else add as normal
        else:
            nameweightdict[namestring] = w

    # convert [[name1, name2, weight],...] to [(index1, index2, weight),...] so it works with add_edge_list()
    iedgelist = []
    for i in nameweightdict.keys():
        n1, n2 = i.split(sep="-", maxsplit=1)  # splits edge name back into constituent node names
        w = nameweightdict[i]
        i1 = nodeindexmap[n1]  #fetch index mapping by name
        i2 = nodeindexmap[n2]  #fetch index mapping by name

        newtuple = (i1, i2, w)
        iedgelist.append(newtuple)

    g.add_edge_list(iedgelist, eprops=[weight])

    # make properties internal before saving internal
    g.vertex_properties["vname"] = vname
    g.edge_properties["weight"] = weight

    g.save(mlname)


# creates hic network from .pairs file
def create_graph_from_pairs(pairsfile, outfile, binsize=250000, verbose=False, outtsv=None):
    # instantiate list of chromosomes
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    if verbose:
        print("getting chromosome sizes...")
    # get chrom sizes from file
    chromsizedict = {}
    with open(pairsfile, "r") as pf:
        for line in pf:
            if line.startswith("#chromsize:"):  # check line for chromsize identifier
                splitline = line.split()
                if splitline[1] in chrlist:  # if chrom is in list we defined above, add that size to chromsizedict
                    chromsizedict[splitline[1]] = int(splitline[2])
                else:  # otherwise skip
                    continue

            elif line.startswith("#columns:"):  # stop file reading once columns are reached
                break
            else:
                continue

    # create empty adjacency matrix according to binsize
    # create {chrompos:index} dict where chrompos is str of format "chromosome:bin#" and index is int index of that node
    startindexdict = {}  # dict for starting starting indices of each chromosome
    nodenamedict = {}  # dict for mapping node# to node name ("chr#:start-stop")
    nodecounter = 0

    # get count of all nodes and nodes for each chromosome
    if verbose:
        print("getting node counts for each chromosome...")
    for chr in chrlist:
         startindexdict[chr] = nodecounter
         nodenamedict[chr] = []
         for i in range(0, math.floor(chromsizedict[chr]/binsize)):  # divide chromosome length by binsize to know bin
            nodecounter = nodecounter + 1  # increment total node counter by one
            # get name of node ("chr#:start-stop") and add to nodenamedict
            numchrombins = math.ceil(chromsizedict[chr] / binsize)
            startpos = i * binsize
            endpos = (i + 1) * binsize - 1
            nodename = chr + ":" + str(startpos) + "-" + str(endpos)
            nodenamedict[chr].append(nodename)
    if verbose:
        print(str(nodecounter-1) + " nodes")

    # create contact matrix
    # create empty ndarray for storing hic contact adjacency matrix
    contactmatrix = np.ndarray((nodecounter, nodecounter))
    if verbose:
        print("populating contact matrix...")
    with open(pairsfile, "r") as pf:  # parse each line for contact information
        for line in pf:
            if line.startswith("#"):  # skip all header lines
                continue
            else:
                chr1 = line.split()[1]
                pos1 = int(line.split()[2])
                chr2 = line.split()[3]
                pos2 = int(line.split()[4])

                try:
                    chr1bin = math.floor(pos1 / binsize) + startindexdict[chr1]  # matrix index for first contact read
                    chr2bin = math.floor(pos2 / binsize) + startindexdict[chr2]  # matrix index for second contact read
                except KeyError as ke:
                    if verbose:
                        print(ke)
                    continue

                contactmatrix[chr1bin, chr2bin] = contactmatrix[chr1bin, chr2bin] + 1
                contactmatrix[chr2bin, chr1bin] = contactmatrix[chr2bin, chr1bin] + 1  # do both ways so symmetrical mat

    if verbose:
        print(contactmatrix.shape)

    if outtsv:
        if verbose:
            print("saving contact matrix tsv...")
        np.savetxt(outtsv, contactmatrix)  # save as tsv before graph creation

    # create graph
    if verbose:
        print("greating graph vertices...")
    # add all vertices from indices and give them names
    g = Graph()
    vname = g.new_vp("string")
    for chr in chrlist:
        for i in range(0, len(nodenamedict[chr])):  # divide chromosome length by binsize to know bin
            v1 = g.add_vertex()
            vname[v1] = nodenamedict[chr][i]  # accesses node name from list of names in nodenamedict

    # add all edges using values from contactmatrix
    if verbose:
        print("creating graph edges...")
    contactcount = g.new_ep("int")
    for i in range(0,nodecounter):
        for j in range(i, nodecounter):  # only traverse the upper triangle of matrix
            if contactmatrix[i, j] > 0:
                v1 = g.vertex(i)
                v2 = g.vertex(j)
                e1 = g.add_edge(v1, v2)
                contactcount[e1] = contactmatrix[i, j]

    # make name and contact graph properties internal
    if verbose:
        print("saving graph file...")
    g.vp["vname"] = vname
    g.ep["raw_contacts"] = contactcount
    g.ep["weight"] = contactcount  # also save as weight property for robustness/compatability
    # save graph
    g.save(outfile)


# uses the fanc api to retrieve TADs
def sam_to_hic(sam1, sam2, outname, outfolder=".", genomefile="hg38.fa", restriction_enzyme="HindIII", nthreads=1):
    # load and sort sam files
    s1 = fanc.load(sam1)
    s2 = fanc.load(sam2)
    sorted_s1 = sort_natural_sam(s1)
    sorted_s2 = sort_natural_sam(s2)

    # load genome regions from genome fasta
    fragments = genome_regions(genomefile, restriction_enzyme=restriction_enzyme)

    # filter reads
    quality_filter = QualityFilter(30, mask='MAPQ')
    uniqueness_filter = UniquenessFilter(strict=True, mask='unique')

    # generate pairs
    pairs_folder = mkdir(os.path.join(outfolder, 'pairs'))
    pairsname = outname + ".pairs"
    pairs = generate_pairs(sorted_s1, sorted_s2, fragments,
                           read_filters=(quality_filter, uniqueness_filter),
                           output_file=os.path.join(pairs_folder, pairsname),
                           check_sorted=True, threads=nthreads)

    # filter pairs
    sl_filter = SelfLigationFilter(mask='self-ligation')
    pairs.filter(sl_filter)

    # obtain hic
    hic_folder = mkdir(os.path.join(outfolder, 'hic'))
    hicname = outname + ".hic"
    hic_file = os.path.join(hic_folder, hicname)
    hic = pairs.to_hic(file_name=hic_file)

    # filter hic
    lc_filter = LowCoverageFilter(binned_hic, rel_cutoff=0.2)  # filter out low coverage contacts
    dg_filter = DiagonalFilter(binned_hic)  # set diagonal to 0
    binned_hic.filter(lc_filter)
    binned_hic.filter(dg_filter)
    binned_hic.run_queued_filters()

    # ICE balance
    ice_balancing(binned_hic)


# uses fanc hic file from above function or other the fanc api function to retrieve TADs
def tads_from_hic():
    print()


def check_contact(binlist, newchrom, newbin):  # listofchromatin contacts from SAMtoGraph,
    contactflag = 0  # flag raised if adjacency is already in list, so it is not appended

    if not binlist:  # check if list of adjacencies is currently empty
        binlist = [[newchrom, newbin, 1]]
    else:
        for i in binlist:  # iterate through the list of adjacencies, +1 weight if found and flag raised
            if newchrom == i[0] and newbin == i[1]:
                i[2] = i[2] + 1
                contactflag = 1  # set flag so a new contact is not appended
                break  # break out of loop once appropriate bin is found
            else:
                continue

        if contactflag:  # if flag raised, do not append contact to list
            pass
        else:
            binlist.append([newchrom, newbin, 1])  # otherwise, do append

    return binlist


# utility function for writing adjlist to file if not being piped directly into hic_adjlist_to_graphml()
def save_adjlist(adjlist, outfile):
    # chromosome list for mapping adjlist 1st index to chromosome name
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    try:
        with open(outfile, "w") as f:
            for chrom in range(len(adjlist)):
                chromname = chrlist[chrom]
                for thisbin in range(len(adjlist[chrom])):
                    for adjacency in adjlist[chrom][thisbin]:
                        writestring = str(chromname) + ":" + str(thisbin) + "\t" + str(adjacency[0]) + ":" + \
                                      str(adjacency[1]) + "\t" + \
                                      str(adjacency[2]) + "\n"  # writes in format: chr1:bin1 \t chr2:bin2 \t weight \n
                        f.write(writestring)  # write line to file
        print("Adjacency list written...")
    except TypeError:
        print("WARNING: TypeError prevented adjlist from writing to file...")


# takes saved adjlist .tsv file and converts it to
def adjlist_file_to_graphml(adjlistfile, outfile="out.graphml"):
    with open(adjlistfile, "r") as f:

        chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
                   'chr12',
                   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
                   'chrY']

        g = Graph(directed=False)

        # create adjacency list to then import into graph with add_edgelist()
        contactlist = []  # for storing contacts
        nodelist = []  # for storing  unique node names, to later be  enumerated
        # initialize graph properties
        vname = g.new_vertex_property("string")  # vp map for node names (chromosome-bin pairs)
        weight = g.new_edge_property("double")  # ep map for edge weight
        # make properties internal
        g.vertex_properties["vname"] = vname
        g.edge_properties["weight"] = weight

        for line in f:
            nodeu, nodev, strweight = line.split()
            floatweight = float(strweight)
            # append nodes to nodelist and contact info to contactlist
            nodelist.append(nodeu)
            nodelist.append(nodev)
            contactlist.append([nodeu, nodev, floatweight])  # is list of lists now, but will later be list of tuples

        # once nodes and edges have been established, enumerate list of nodes and create graph using (index1,index2,weight)
        nodelist = list(set(nodelist))  # keep only unique entries in list of nodenames
        nlist = list(enumerate(nodelist))  # enumerate list for passing of index to graph-tool with add_edge_list()
        nodeindexmap = {}  # empty dict for mapping value:index pairs
        g.add_vertex(n=len(nodelist))  # add all vertices to graph first

        for n in nlist:
            nodeindexmap[n[1]] = n[0]  # add node index:name mapping
            vname[n[0]] = n[1]  # add name to vertex properties

        # prune edge list so a->b, b->a duplicates are combined and their weights are summed
        # uses a dictionary for ease of ooking up/ modifying weight value
        nameweightdict = {}
        for i in contactlist:
            n1, n2, w = i  # unpack "tuple" (actually a list rn but)
            # initialize strings for forward and reverse edge direction
            namestring = ""
            namestring = namestring.join([n1, "-", n2])
            rvrstring = ""
            rvrstring = rvrstring.join([n2, "-", n1])  # reverse string in case some edges are double encoded

            # if one name is already in the dictionary, update weight instead of adding
            if namestring in nameweightdict or rvrstring in nameweightdict:
                if namestring in nameweightdict:
                    nameweightdict[namestring] = nameweightdict[namestring] + w
                if rvrstring in nameweightdict:
                    nameweightdict[rvrstring] = nameweightdict[rvrstring] + w
            # else add as normal
            else:
                nameweightdict[namestring] = w

        # convert [[name1, name2, weight],...] to [(index1, index2, weight),...] so it works with add_edge_list()
        iedgelist = []
        for i in nameweightdict.keys():
            n1, n2 = i.split(sep="-", maxsplit=1)  # splits edge name back into constituent node names
            w = nameweightdict[i]
            i1 = nodeindexmap[n1]  # fetch index mapping by name
            i2 = nodeindexmap[n2]  # fetch index mapping by name

            newtuple = (i1, i2, w)
            iedgelist.append(newtuple)

        g.add_edge_list(iedgelist, eprops=[weight])

    g.save(outfile)


# creates an edgelist from a .hic file to later be converted into a graph format (saves memory)
def edgelist_from_hic(hicfile, outfile, binsize=100000, norm="None", contact_type="oe"):

    hic = hicstraw.HiCFile(hicfile)
    chrom_list = []

    for chrom in hic.getChromosomes():
      chrom_list.append(chrom.name)

    matrix_object = hic.getMatrixZoomData('chr1', 'chr2', "observed", "None", "BP", 250000)  # alt way of getting contacts
    print(matrix_object)

    with open(outfile, "w") as f:
        for first_chrom in tqdm(chrom_list, total=len(chrom_list)):
            if first_chrom != "All":
                for second_chrom in chrom_list:
                    if second_chrom != "All":
                        result = hicstraw.straw(contact_type, norm, hicfile, first_chrom, second_chrom, 'BP', binsize)
                        for i in range(len(result)):
                            # v1 = first_chrom + ":" + str(int(int(result[i].binX) / binsize))
                            # v2 = second_chrom + ":" + str(int(int(result[i].binY) / binsize))
                            v1 = first_chrom + ":" + str(int(int(result[i].binX)))
                            v2 = second_chrom + ":" + str(int(int(result[i].binY)))
                            f.write("{0}\t{1}\t{2}\n".format(v1, v2, result[i].counts))


# reads edgelist .tsv file from edgelist_from_hic() and converts it into a .gt graph file
def graph_from_hic_edgelist(edgelist_file, outfile, interchrom_only=False):

    g = Graph(directed=False)
    weight = g.new_edge_property("double")  # ep map for edge weight
    g.edge_properties["weight"] = weight  # make weight an internal property
    edge_list = []
    vname_index_dict = {}
    index_counter = 0

    # count lines in file
    with open(edgelist_file, "rb") as f:
        num_lines = sum(1 for _ in f)

    with open(edgelist_file, "r") as f:
        for line in tqdm(f, total=num_lines):
            v1str, v2str, wstr = line.split("\t")
            if interchrom_only:
                if v1str.split(":")[0] == v2str.split(":")[0]:
                    continue
            if v1str not in vname_index_dict.keys():
                vname_index_dict[v1str] = index_counter
                index_counter = index_counter + 1
            if v2str not in vname_index_dict.keys():
                vname_index_dict[v2str] = index_counter
                index_counter = index_counter + 1
            v1 = vname_index_dict[v1str]
            v2 = vname_index_dict[v2str]

            e = g.add_edge(v1, v2)
            g.ep.weight[e] = float(wstr)

            # edge_list.append((v1, v2, float(wstr)))
        # g.add_edge_list(edge_list, eprops=weight)

        print("Hoi!")
        g.save(outfile)
        print("Hoi!")
    # label all vertices with genome bin ids
    vname = g.new_vertex_property("string")  # vp map for node names (chromosome-bin pairs)
    for i in vname_index_dict.keys():
        vname[vname_index_dict[i]] = i

    g.vertex_properties["vname"] = vname
    print("Hoi!")
    g.save(outfile)
    print("Hoi!")


# function for retrieving all keys with a given value
def get_keys_by_value(adict, value):
    keylist = []
    itemlist = adict.items()
    for item in itemlist:
        for id in item[1]:
            if id == value:
                keylist.append(item[0])
    return keylist


# annotates a given hic network with genes and go terms found in gencodefile and mapfile respectively
def genes_go_annotate(graphfile, mapfile="HUMAN_9606_idmapping_selected.tsv", gencodefile="gencode_pcgenes.csv",
                          binsize=250000, outfile="annotatedgonet.graphml", go_obo="go-basic.obo", binfile=None):
    start_time = time.time()

    go = obo_parser.GODag(go_obo)  # initialize go file

    gencode = pd.read_csv(gencodefile)  # load gencode table
    mapping = pd.read_csv(mapfile, sep='\t')  # load mapping table

    # create attribute dictionaries, dict name will be attribute name in graph
    goterms = {}  # key is node name (chr[#]:[bin]), value is list of goterms
    genes = {}  # key is node name (chr[#]:[bin]), value is list of genes

    try:
        if binfile is None:
            # create new vname prop to match old style (chr#:chromnodenumber vs chr#:startpos-endpos)
            g = load_graph(graphfile)
            newvname = g.new_vertex_property("string")
            for v in g.vertices():
                vn = g.vp.vname[v]
                chrom = vn.split(":")[0]
                startpos = int(vn.split(":")[1].split("-")[0])
                newvname[v] = chrom + ":" + str(int(startpos/binsize))
                print(newvname[v])

            # for each gene, take centroid between start and end pos, then convert pos to bin. annotate that bin
            print("creating gene annotations...")
            print("--- %s seconds since start ---" % (time.time() - start_time))
            for index, row in gencode.iterrows():
                chrom = row['seqid']
                startpos = int(row['start'])
                endpos = int(row['end'])
                centroid = int(startpos + endpos / 2)
                centrebin = int(centroid / binsize)

                # for each node to be annotated, check if key exists, then create/append new gene as necessary
                nodename = chrom + ':' + str(centrebin)
                if nodename in genes:  # if key already exists, append new gene
                    genes[nodename].append(row['id'])
                else:  # else make new gene list that contains gene
                    genes[nodename] = [row['id']]

        else:  # if binfile is included, use that (should be bed format chr#\wstart\wend)
            # for each bin, take centroid between start and end pos, then convert pos to bin. annotate that bin
            print("creating gene annotations...")
            print("--- %s seconds since start ---" % (time.time() - start_time))

            # make a dict for bin ranges
            bindict = {"chr1": [], "chr2": [], "chr3": [], "chr4": [], "chr5": [], "chr6": [], "chr7": [], "chr8": [],
                       "chr9": [], "chr10": [], "chr11": [], "chr12": [], "chr13": [], "chr14": [], "chr15": [],
                       "chr16": [],
                       "chr17": [], "chr18": [], "chr19": [], "chr20": [], "chr21": [], "chr22": [], "chrX": [],
                       "chrY": [],
                       "chrM": []}
            # populate it with bin start and end positions
            f = open(binfile, "r")
            for line in f:
                chrom, startpos, endpos = line.split()
                binbounds = str(startpos) + "-" + str(endpos)
                bindict[chrom].append(binbounds)
            f.close()

            # iterate through all protein coding genes and add them to appropriate bins
            for index, row in gencode.iterrows():
                chrom = row['seqid']
                startpos = int(row['start'])
                endpos = int(row['end'])
                centroid = int(startpos + endpos / 2)
                genename = row['id']

                indexcounter = 0  # counts iters to keep track of bin index, which is used in vname ("[chr]:[bindex]")
                for whichbin in bindict[chrom]:
                    binstart, binend = whichbin.split("-")  # get bin start and end positions from name
                    if int(binstart) <= centroid <= int(binend):  # if centroid is in bin
                        nodename = chrom + ':' + str(indexcounter)  # nodename set to match vname
                        if nodename in genes:  # if key already exists, append new gene
                            genes[nodename] = genes[nodename] + [genename]
                        else:  # else make new gene list that contains gene
                            genes[nodename] = [genename]
                    indexcounter = indexcounter + 1  # increments index

        # for each node, for each gene, add associated GO terms to new dictionary with node as key, GO term as value
        print("creating GO annotations...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        for ind, ro in tqdm(mapping.iterrows()):  # iterate through mapping table
            ids = str(ro.Ensembl)  # split multiple ids into individual ones
            ids = ids.split("; ")
            for emblid in ids:  # iterate through those ensembl ids
                keys = get_keys_by_value(genes, emblid)  # get all nodes with that gene
                for i in keys:  # iterate through that
                    termstring = str(ro.GO)
                    terms = termstring.split("; ")  # coerce string (list of go terms) to list
                    if i not in goterms:  # if no go terms yet, initialize list
                        goterms[i] = terms
                    else:  # else append to list
                        goterms[i] = goterms[i] + terms

        # for each node, for each GO term add all parent terms to dictionary
        print("adding parent GO terms...")
        print("--- %s seconds since start ---" % (time.time() - start_time))

        for node in tqdm(goterms.keys()):
            try:
                ogterms = goterms[node]  # original list of terms
                if 'nan' in ogterms:
                    ogterms.remove('nan')  # remove NaNs
                allterms = ogterms  # list for all terms to be added to, starting with oglist
                for term in ogterms:
                    if term != 'nan':  # skip NaNs
                        rec = go[term]
                        parents = list(rec.get_all_parents())
                        allterms = allterms + parents
                allterms = list(set(allterms))  # removes duplicates from list
                goterms[node] = allterms
            except KeyError as keyerror:
                print("keyerror: ")
                print(keyerror)  # prints GOTERMS not found in file

            except TypeError:  # ignore nans if they still exist
                pass

        print("adding annotations to graph...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        # initialize new vertex property objects
        geneprop = g.new_vertex_property("vector<string>")
        goprop = g.new_vertex_property("vector<string>")

        if binfile is None:
            x=0  # TODO remove
            # loop iterates through all gene annotations and adds them to geneprop
            for k in genes.keys():
                try:
                    v = g.vertex_index[find_vertex(g, newvname, str(k))[0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                    geneprop[v] = genes[k]
                except IndexError:
                    x = x+1  # REMOVE
                    # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                    print("no available vertex for " + str(k) + ". skipping...")
                    print(newvname[g.vertex(x)])  # REMOVE

            # loop iterates through all go annotations and adds them to goprop
            for k in goterms.keys():
                try:
                    v = g.vertex_index[find_vertex(g, newvname, str(k))[
                        0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                    goprop[v] = goterms[k]
                except IndexError:
                    # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                    print("no available vertex for " + str(k) + ". skipping...")

        else:
            # loop iterates through all gene annotations and adds them to geneprop
            for k in genes.keys():
                try:
                    v = g.vertex_index[find_vertex(g, g.vp["vname"], str(k))[
                        0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                    geneprop[v] = genes[k]
                except IndexError:
                    # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                    print("no available vertex for " + str(k) + ". skipping...")

            # loop iterates through all go annotations and adds them to goprop
            for k in goterms.keys():
                try:
                    v = g.vertex_index[find_vertex(g, g.vp["vname"], str(k))[
                        0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                    goprop[v] = goterms[k]
                except IndexError:
                    # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                    print("no available vertex for " + str(k) + ". skipping...")

        # make property maps internal
        g.vertex_properties["genes"] = geneprop  # networkx naming convention kept for back compatibility
        g.vertex_properties["goterms"] = goprop  # networkx naming convention kept for back compatibility

        print("writing graph to file...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        g.save(outfile)

    except MemoryError as memerror:
        print(memerror)

    print("Annotated graph saved as " + outfile)
    print("--- %s seconds since start ---" % (time.time() - start_time))


# annotates a given hic network with genes and go terms found in gencodefile and mapfile respectively
# DEPRECATED due to new way of using vname for node position
def genes_go_annotate_old(graphfile, mapfile="HUMAN_9606_idmapping_selected.tsv", gencodefile="gencode_pcgenes.csv",
                      m=40000, outfile="annotatedgonet.graphml", go_obo="go-basic.obo", binfile=None):

    start_time = time.time()

    go = obo_parser.GODag(go_obo)  # initialize go file

    gencode = pd.read_csv(gencodefile)  # load gencode table
    mapping = pd.read_csv(mapfile, sep='\t')  # load mapping table

    # create attribute dictionaries, dict name will be attribute name in graph
    goterms = {}  # key is node name (chr[#]:[bin]), value is list of goterms
    genes = {}  # key is node name (chr[#]:[bin]), value is list of genes

    try:
        if binfile is None:
            # for each bin, take centroid between start and end pos, then convert pos to bin. annotate that bin
            print("creating gene annotations...")
            print("--- %s seconds since start ---" % (time.time() - start_time))
            for index, row in gencode.iterrows():
                chrom = row['seqid']
                startpos = int(row['start'])
                endpos = int(row['end'])
                centroid = int(startpos + endpos / 2)
                centrebin = int(centroid / m)

                # for each node to be annotated, check if key exists, then create/append new gene as necessary
                nodename = chrom + ':' + str(centrebin)
                if nodename in genes:  # if key already exists, append new gene
                    genes[nodename].append(row['id'])
                else:  # else make new gene list that contains gene
                    genes[nodename] = [row['id']]

        else:  # if binfile is included, use that (should be bed format chr#\wstart\wend)
            # for each bin, take centroid between start and end pos, then convert pos to bin. annotate that bin
            print("creating gene annotations...")
            print("--- %s seconds since start ---" % (time.time() - start_time))

            # make a dict for bin ranges
            bindict = {"chr1": [], "chr2": [], "chr3": [], "chr4": [], "chr5": [], "chr6": [], "chr7": [], "chr8": [],
                       "chr9": [], "chr10": [], "chr11": [], "chr12": [], "chr13": [], "chr14": [], "chr15": [],
                       "chr16": [],
                       "chr17": [], "chr18": [], "chr19": [], "chr20": [], "chr21": [], "chr22": [], "chrX": [],
                       "chrY": [],
                       "chrM": []}
            # populate it with bin start and end positions
            f = open(binfile, "r")
            for line in f:
                chrom, startpos, endpos = line.split()
                binbounds = str(startpos) + "-" + str(endpos)
                bindict[chrom].append(binbounds)
            f.close()

            # iterate through all protein coding genes and add them to appropriate bins
            for index, row in gencode.iterrows():
                chrom = row['seqid']
                startpos = int(row['start'])
                endpos = int(row['end'])
                centroid = int(startpos + endpos / 2)
                genename = row['id']

                indexcounter = 0  # counts iters to keep track of bin index, which is used in vname ("[chr]:[bindex]")
                for whichbin in bindict[chrom]:
                    binstart, binend = whichbin.split("-")  # get bin start and end positions from name
                    if int(binstart) <= centroid <= int(binend):  # if centroid is in bin
                        nodename = chrom + ':' + str(indexcounter)  # nodename set to match vname
                        if nodename in genes:  # if key already exists, append new gene
                            genes[nodename] = genes[nodename] + [genename]
                        else:  # else make new gene list that contains gene
                            genes[nodename] = [genename]
                    indexcounter = indexcounter + 1  # increments index

        # for each node, for each gene, add associated GO terms to new dictionary with node as key, GO term as value
        print("creating GO annotations...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        for ind, ro in tqdm(mapping.iterrows()):  # iterate through mapping table
            ids = str(ro.Ensembl)  # split multiple ids into individual ones
            ids = ids.split("; ")
            for emblid in ids:  # iterate through those ensembl ids
                keys = get_keys_by_value(genes, emblid)  # get all nodes with that gene
                for i in keys:  # iterate through that
                    termstring = str(ro.GO)
                    terms = termstring.split("; ")  # coerce string (list of go terms) to list
                    if i not in goterms:  # if no go terms yet, initialize list
                        goterms[i] = terms
                    else:  # else append to list
                        goterms[i] = goterms[i] + terms

        # for each node, for each GO term add all parent terms to dictionary
        print("adding parent GO terms...")
        print("--- %s seconds since start ---" % (time.time() - start_time))

        for node in tqdm(goterms.keys()):
            try:
                ogterms = goterms[node]  # original list of terms
                if 'nan' in ogterms:
                    ogterms.remove('nan')  # remove NaNs
                allterms = ogterms  # list for all terms to be added to, starting with oglist
                for term in ogterms:
                    if term != 'nan':  # skip NaNs
                        rec = go[term]
                        parents = list(rec.get_all_parents())
                        allterms = allterms + parents
                allterms = list(set(allterms))  # removes duplicates from list
                goterms[node] = allterms
            except KeyError as keyerror:
                print("keyerror: ")
                print(keyerror)  # prints GOTERMS not found in file

            except TypeError:  # ignore nans if they still exist
                pass

        print("adding annotations to graph...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        g = load_graph(graphfile)
        # initialize new vertex property objects
        geneprop = g.new_vertex_property("vector<string>")
        goprop = g.new_vertex_property("vector<string>")

        # loop iterates through all gene annotations and adds them to geneprop
        for k in genes.keys():
            try:
                v = g.vertex_index[find_vertex(g, g.vp["vname"], str(k))[0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                geneprop[v] = genes[k]
            except IndexError:
                # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                print("no available vertex for " + str(k) + ". skipping...")

        # loop iterates through all go annotations and addds them to goprop
        for k in goterms.keys():
            try:
                v = g.vertex_index[find_vertex(g, g.vp["vname"], str(k))[0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                goprop[v] = goterms[k]
            except IndexError:
                # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                print("no available vertex for " + str(k) + ". skipping...")

        # make property maps internal
        g.vertex_properties["genes"] = geneprop  # networkx naming convention kept for back compatibility
        g.vertex_properties["goterms"] = goprop  # networkx naming convention kept for back compatibility

        print("writing graph to file...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        g.save(outfile)

    except MemoryError as memerror:
        print(memerror)

    print("Annotated graph saved as " + outfile)
    print("--- %s seconds since start ---" % (time.time() - start_time))


# use mapfile and genes annotated to nodes to annotate nodes with REACTOME pathways
def reactome_annotate(graph, mapfile="Ensembl2Reactome_All_Levels.txt", outfile=None):
    listoflists = []  # initialize list for coersion to df

    # read lines and append to listoflists
    with open(mapfile, "r") as f:

        for line in f.readlines():
            splitline = line.split(sep="\t")
            splitline[-1] = splitline[-1][:-1]  # remove endline character
            listoflists.append(splitline)

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    # filter table so only human genes included
    df = pd.DataFrame.from_records(listoflists)  # coerce to dataframe
    df = df.rename({0: "EMBLid", 1: "Pathway", 2: "URL", 3: "Notes", 4: "Source", 5: "Species"}, axis=1)
    df = df.loc[df["Species"] == "Homo sapiens"]
    # print(df.loc[df["EMBLid"] == "ENSG00000129038"])
    # print(find_vertex(g, g.vp.genes, 'ENSG00000129038'))

    reactomepathways = g.new_vp("vector<string>")

    # iterate over vertices, extract gene list, use genes to annotate with pathways
    for v in g.vertices():
        genelist = g.vp.genes[v]

        # create list of all reactome paths associated with all genes on a node
        reactomelist = []
        for gene in genelist:
            rows = df.loc[df["EMBLid"] == gene]
            for index, row in rows.iterrows():
                reactomelist.append(row["Pathway"])

        # assign full list to node as vertex property
        reactomepathways[v] = reactomelist

    if outfile == None:
        if type(graph) == str:
            filename = graph[:-3] + "_REACTOMEannotated.gt"
            g.save(filename)
        else:
            filename = "REACTOMEannotatedGraph.gt"
            g.save(filename)
    else:
        g.save(outfile)


# converts gencode gff3 file to parsed, filtered table
def parse_gff3(infile, outfile='./gencode_pcgenes.csv'):
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    mylist = [[]]  # empty list for later coersion to data frame

    with open(infile) as myfile:
        lines = myfile.readlines()
    for line in lines:  # goes line by line and coerces final table to data frame
        if not line.startswith('#'):  # ignores header line
            mylist.append(line.split()[0:9])
    mycolumns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gencode = pd.DataFrame(mylist, columns=mycolumns)

    gencode = gencode.where(gencode.type == "gene")  # only include genes in table
    gencode = gencode.where(gencode.seqid.isin(chrlist))  # only include nuclear chromosomes
    gencode['id'] = gencode.attributes.str.split(';').str.get(0)
    gencode['gene_type'] = gencode.attributes.str.split(';').str.get(2)
    gencode = gencode.dropna()
    gencode = gencode.where(gencode.gene_type == "gene_type=protein_coding")  # only include protein coding genes
    gencode = gencode.dropna()
    gencode['id'] = gencode.id.str.split('=').str.get(1)  # remove 'gene_id=' from id
    gencode['id'] = gencode.id.str.split('.').str.get(0)  # remove decimals from end of gene id

    gencode.to_csv(outfile)



# first though, this function makes the COO from a graph
def make_graph_coo(graph, outfile="graphCOO.tsv"):

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    f = open(outfile, "w")  # TODO name this properly? or delete after running?

    # get edgelist
    # convert to COO and store in file
    for s, t, w in g.iter_edges([g.ep.weight]):  # source, target, weight
        f.write(str(s) + " " + str(t) + " " + str(w) + "\n")
    f.close()


def read_edge_list(file_path):
    """
    Reads an edge list from a TSV file and converts it into a sparse matrix.

    Parameters:
    - file_path (str): Path to the TSV file containing the edge list.

    Returns:
    - scipy.sparse.coo_matrix: Sparse matrix representation of the edge list.
    """
    df = pd.read_csv(file_path, sep='\t')
    sources = df['source']
    targets = df['target']
    weights = df['weight']

    # Convert sources and targets to categorical types for consistent indexing
    sources = pd.Categorical(sources)
    targets = pd.Categorical(targets)

    # Create a sparse matrix
    matrix = coo_matrix((weights, (sources.codes, targets.codes)))

    return matrix


# filters out spurious contacts before normalization
def prenorm_filter(graph, outfile, min_contacts=2):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    edges_to_remove = []
    for e in g.edges():
        if g.ep.weight[e] < min_contacts:
            edges_to_remove.append(e)

    for e in edges_to_remove:
        g.remove_edge(e)

    g.save(outfile)


# extracts the largest connected component from the graph and saves it as a new file
def save_largest_component_graph(graph, outfile="graphLCC.gt"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    fg = extract_largest_component(g, prune=True)

    fg.save(outfile)


# performs ICE normalization (imakaev et al.) on COO matrix
def ice_balance(infile, outfile):
    # Load data from file
    data = pd.read_csv(infile, sep=' ', header=None, names=['source', 'target', 'weight'],
                       dtype = {'source': np.int64, 'target': np.int64, 'weight': np.float128})

    # Create sparse matrix, concatenating with reversed ro/col to ensure symmetrical
    row = np.concatenate((data['source'].values, data['target'].values))
    col = np.concatenate((data['target'].values, data['source'].values))
    values = np.concatenate((data['weight'].values, data['weight'].values))
    n = max(np.max(row), np.max(col)) + 1  # Assuming 0-based index, adjust accordingly
    mat = coo_matrix((values, (row, col)), shape=(n, n))

    # do ice
    normed_mat = normalization.ICE_normalization(mat).toarray()
    del mat
    print("normalization done. Writing to file...")

    tril_normed_mat = np.tril(normed_mat, k=-1)  # get triangular matrix without the main diagonal
    del normed_mat

    writelist = []
    sources, targets = tril_normed_mat.nonzero()
    for src, tgt in tqdm(zip(sources.tolist(), targets.tolist()), total=len(sources)):
        writelist.append(str(src) + " " + str(tgt) + " " + str(tril_normed_mat[src, tgt]) + "\n")

    with open(outfile, "w") as f:
        f.writelines(writelist)


# updates edgeweights from a COO file. Mostly for use with ice_balance() or other normalization methods
# actually creates a new graph rather than updating old one, as updating turns out to be very slow
def update_weights_from_coo(graph, coo, outgraph):
    # update edge weights with normalized values
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    # create new graph with vertices from original graph
    g2 = Graph(directed=False)
    g2vname = g2.new_vertex_property("string")
    for v in g.vertices():
        v1 = g2.add_vertex()
        g2vname[v1] = g.vp.vname[v]  # copies vname (str chromosome position) to copied vertex in new graph
    g2.vp["vname"] = g2vname
    print("Vertices added to new graph")

    eweight_dict = g2.new_edge_property("double")
    # draw edges based on new weights
    with open(coo, "r") as cf:
        lc = len(cf.readlines())
    with open(coo, "r") as cf:
        for line in tqdm(cf, total=lc):
            splitline = line.split()
            src = int(splitline[0])
            tgt = int(splitline[1])
            val = float(splitline[2])

            e1 = g2.add_edge(g2.vertex(src), g2.vertex(tgt))
            eweight_dict[e1] = val

    g2.ep["weight"] = eweight_dict  # make weight property internal

    # copy mods graph property
    try:
        g2.gp["mods"] = g.gp.mods
    except ValueError:
        print("Could not add mods graph property to updated graph")
    # write new graph file
    g2.save(outgraph)


def remove_self_edges(graph, outgraph):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    vcount = count_vertices(g)
    g.set_fast_edge_removal()

    # create list of edge tuples for all adges with identical source and target vertices
    sedgelist = [(x, x) for x in np.arange(0,vcount)]

    for s,t in sedgelist:
        sedge = g.edge(g.vertex(s),g.vertex(t))
        if not isinstance(sedge, type(None)):
            g.remove_edge(sedge)

    # add diagonal removal to list of mods
    try:
        g.gp.mods.append("Self-edges removed")
    except AttributeError:
        mods = g.new_gp("vector<string>")
        g.gp["mods"] = mods
        g.gp.mods.append("Self-edges removed")

    g.save(outgraph)


# make a dictionary of 'GOterm': [nodes with that term]
def make_go_dict(graph, prop="goterms", writefile=None):
    godict = {}  # initialize dict
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    for v in g.vertices():  # iterate through nodes
        if prop == "goterms":
            terms = g.vp.goterms[v]
        elif prop == "reactomepathways":
            terms = g.vp.reactomepathways[v]

        try:
            for term in terms:  # iterate through terms
                if term in godict:
                    godict[term].append(g.vertex_index[v])  # append new node to list if it already exists
                else:
                    godict[term] = [g.vertex_index[v]]  # else make a new list for terms
        except KeyError as keyerror:
            pass
    print(prop + " dictionary created")

    if writefile is not None:  # write to named file if requested
        with open(writefile, "w") as f:
            for k in godict.keys():
                writestring = str(k) + " " + " ".join([str(x) for x in godict[k]]) + "\n"
                f.write(writestring)

    return godict


# writes go dict to file from annotated graph
def write_go_dict(graph, outfile):
    # parse graph input
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    #write to file
    with open(outfile, "w") as f:
        godict = make_go_dict(g)

        for k in godict.keys():
            writestring = str(k) + " " + " ".join([str(x) for x in godict[k]]) + "\n"
            f.write(writestring)


# make a histogram where y is number of GOTERMs found on x number of nodes
def plot_goterms_per_node_histogram(adict, outfile="GOhistogram.pdf", xlim=250, ylim=250):
    mydict = {}
    dictogram = np.array(list(adict.items()))

    for row in tqdm(dictogram):
        if len(row[1]) in mydict:
            mydict[len(row[1])] += 1
        else:
            mydict[len(row[1])] = 1

    f = plt.bar(list(mydict.keys()), mydict.values())
    plt.xlim(0, xlim)
    plt.ylim(0, ylim)
    plt.show()
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()

    for item, value in mydict.items():
        print(item, ":", value)


# update all weights in graphs so that lower weight = closer proximity so dijkstra does what we want
def invert_weights(graphfile, outfile="inverted.graphml"):
    g = load_graph(graphfile)
    newweights = g.new_edge_property("double")  # new vp for inverted weights

    for s, t, w in g.iter_edges([g.ep.weight]):  # source, target, weight:
        weight = w
        if 0 > weight > 1:
            print("Please make sure the graph has been normalized before invoking invert_weights()")
            sys.exit()
        else:
            weight = 1 - weight
            newweights[g.edge(s, t)] = weight

    g.ep["weight"] = newweights  # overwrite weight vp with new weights and make internal before saving

    g.save(outfile)


# for ICE normalization, scale then invert the weights
def invert_scale_weights(graphfile, wannotation="weight", outfile="invertScaled.graphml"):
    g = load_graph(graphfile)

    if wannotation in ["weight", "raw_weight"]:
        newweights = g.new_edge_property("double")  # initialize new property map for updated edge weights
        oldweights = g.ep[wannotation]  # store old weight property map for posterity

        # Scale weights to a range [0, 1]
        if np.any(0 > oldweights.a) or np.any(oldweights.a > 1):  # only apply scaling if weights are not already in range
            owmax = np.max(oldweights.a)
            owmin = np.min(oldweights.a)
            scaledat = (oldweights.a - owmin) / (owmax - owmin)

            for e, sw in zip(g.edges(), scaledat):
                newweights[e] = 1 - sw

            g.ep["weight"] = newweights  # overwrite weight vp with new weights and make internal before saving
            g.ep["raw_weight"] = oldweights

            g.save(outfile)

        else:
            print("already scaled between 0 and 1")

    else:
        print("bad weight annotation")


# just do 1/x for ICE normalized values?
def simple_invert_weights(graphfile, wannotation="weight", outfile="modified.graphml", eps=0.0000001, addone=True):
    g = load_graph(graphfile)

    if wannotation in ["weight", "raw_weight"]:
        newweights = g.new_edge_property("double")  # initialize new property map for updated edge weights
        oldweights = g.ep[wannotation]  # store old weight property map for posterity

        if addone:  # add 1 to ensure there are no 0-1 weights before inverting
            for e, ow in zip(g.edges(), oldweights.a):
                ow = ow + 1
                newweights[e] = 1/ow

        else:  # if option not selected still add some small epsilon to avoid divide by zero errors
            for e, ow in zip(g.edges(), oldweights.a):
                if ow != 0:
                    newweights[e] = 1/ow
                else:
                    ow = ow + eps
                    newweights[e] = 1 / ow

        g.ep["weight"] = newweights  # overwrite weight vp with new weights and make internal before saving
        g.ep["raw_weight"] = oldweights

    else:
        print("Invalid weight annotation. You must choose either \"weight\" or \"raw_weight\".")

    # add modification to mods list graph property
    try:
        g.gp.mods.append("Annotated with genes")
        g.gp.mods.append("Annotated with goterms")
    except AttributeError:
        mods = g.new_gp("vector<string>")
        g.gp["mods"] = mods
        g.gp.mods.append("Annotated with genes")
        g.gp.mods.append("Annotated with goterms")
    g.save(outfile)


# update all weights in graphs so that lower weight = closer proximity so dijkstra does what we want
# old version for graphs where the weight was previously inverted
def modify_weights(graphfile, wannotation="weight", outfile="modified.graphml"):
    g = load_graph(graphfile)

    if wannotation in ["weight", "raw_weight"]:
        newweights = g.new_edge_property("double")  # initialize new property map for updated edge weights
        oldweights = g.ep[wannotation]  # store old weight property map for posterity

        # Normalize weights to a range [0, 1]
        maxdenom = max(1/(1-oldweights.a))
        print(maxdenom)
        wmax = max(oldweights.a)  # for skipping this divide by zero error during transformation
        print(wmax)

        for e, ow in zip(g.edges(), oldweights.a):
            if ow == wmax:
                newweights[e] = float('inf')
            else:
                newweights[e] = 1/-np.log10((1/(1-ow))/maxdenom)

        # replace infinite values with max weight
        inf_indices = np.where(np.isinf(newweights.a))[0]
        inf_edges = [edge for i, edge in enumerate(g.edges()) if i in inf_indices]
        for e in inf_edges:
            newweights[e] = max(newweights.a[np.where(np.isfinite(newweights.a))])

        # remove after troubleshooting
        print("new weights stats post inf removal")
        print("max: " + str(max(newweights.a)))
        print("med: " + str(np.median(newweights.a)))
        print("min: " + str(min(newweights.a)))
        print("has inf? " + str(np.any(np.isinf(newweights.a))))
        print("has nan? " + str(np.any(np.isnan(newweights.a))))

        g.ep["weight"] = newweights  # overwrite weight vp with new weights and make internal before saving
        g.ep["raw_weight"] = oldweights

    else:
        print("Invalid weight annotation. You must choose either \"weight\" or \"raw_weight\".")

    g.save(outfile)


# converts graphml to gt for faster opening in graph tool
def graphml_to_gt(graph, outfile="graph.gt"):
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
    g.save(outfile, fmt="gt")

# converts gt to graphml for utility
def gt_to_graphml(graph, outfile="graph.graphml"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    g.save(outfile, fmt="graphml")


# splits list of goterms (which is a string) to property map for use with find_vertex()
# Is this necessary? Update: No, but it's a good reference for how property maps work in graph-tool
def map_goterms(graph):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    goterm = g.new_vp("vector<string>")
    for v in g.vertices():
        gostring = g.vp.goterms[v]
        # clean up string before splitting into list
        gostring = gostring.replace("\'", "")
        gostring = gostring.replace(",", "")
        gostring = gostring.replace("[", "")
        gostring = gostring.replace("]", "")
        golist = gostring.split()

        goterm[v] = golist

    if type(graph) == str:
        g.save(graph, fmt="gt")
    else:
        return g


# function that obtains sampling vector from network
# each vertex is added a number of times equal to the number of specified annotations it has
# possible annotations are "goterms", "genes", "reactomepathways"
# TODO dont include goterms above threshold that will be counted
def get_sampling_vector(graph, vmax=None, prop="goterms"):
    if prop == "goterms":
        nndict = make_numnodes_dict(graph)
    elif prop == "reactomepathways":
        nndict = make_numnodes_dict(graph, prop="reactomepathways")

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    samplevec = []
    for v in g.vertices():
        golist = g.vp.goterms[v]

        if vmax:  # if vmax is set, don't append node when annotated with goterm with higher number of nodes
            for i in golist:
                if nndict[i] <= vmax:
                    samplevec.append(g.vertex_index[v])

        else:  # otherwise append node for every term
            for i in golist:
                samplevec.append(g.vertex_index[v])

    return samplevec


# pull a weighted sample from the population
# works by pulling samples randomly from the list, coercing to set to remove duplicates, and repeating until full
def get_random_tpd(distarray, samplevec, vlength):
    samplepop = random.sample(samplevec, vlength)  # get a random sample from our weighted vector

    # until there are vlength unique vertices, keep drawing from weighted vector
    while len(set(samplepop)) != vlength:
        samplepop.extend(random.sample(samplevec, vlength - len(set(samplepop))))
    samplepop = list(set(samplepop))  # condense down to appropriately sized list of unique vals

    combs = combinations(samplepop, 2)  # get all possible combinations for sampled nodes
    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory
    row = list(row)
    col = list(col)
    tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
    return tpd


# top percent pairwise distance
# pull a weighted sample from the population
# works by pulling samples randomly from the list, coercing to set to remove duplicates, and repeating until full
def get_random_tptpd(distarray, samplevec, vlength, tpthresh=0.2):
    samplepop = random.sample(samplevec, vlength)  # get a random sample from our weighted vector

    # until there are vlength unique vertices, keep drawing from weighted vector
    while len(set(samplepop)) != vlength:
        samplepop.extend(random.sample(samplevec, vlength - len(set(samplepop))))

    # TODO assess clustering of nodes and cut out bottom% nodes
    # for speed, put all combs into numpy array, vectorize dist lookup?
    # np.fsum(where=) to get separate arrays for each node, then get sum(tpd) for each before sorting and cutting
    samplepop = list(set(samplepop))  # condense down to appropriately sized list of unique vals
    combs = combinations(samplepop, 2)  # get all possible combinations for sampled nodes
    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory

    # get array of combinations with respective shortest path distances
    distlist = distarray[row, col]
    tosumarray = np.array([[[row, col]],
                           [[distlist, distlist]]])  # create a 2x2xncombs np array with axes source x target x dist
    del distlist

    # for each node, get it's individual tpd and add it to a list to be sorted and trimmed
    tpdlist = []
    for i in range(len(samplepop)):
        mask = tosumarray[0, :, :] == i
        tpd = np.sum(tosumarray[1, :, :], where=mask)  # vectorized addition of PD at all combs for TPD calculation
        tpdlist.append(tpd)

    # add nodes with tpds into data frame, sort it, calculate regular tpd for top% of nodes
    tosortframe = pd.DataFrame({"node": samplepop, "tpd": tpdlist})
    sortedframe = tosortframe.sort_values('tpd', ascending=False)

    numtop = int(vlength * tpthresh)
    if numtop <= 3:  # make sure there are always at least 3 nodes in set
        numtop = 3
        # TODO make this so it can handle 2?

    newnodeslist = sortedframe['node'].tolist()
    newnodeslist = newnodeslist[0:numtop]  # cut down list to only top nodes

    #  recalculate tpd for only top nodes
    combs = combinations(newnodeslist, 2)  # get all possible combinations for sampled nodes
    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory
    row = list(row)
    col = list(col)
    tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation

    return tpd


# generate a vector containing the TPD distribution for a given number of vertices (vlength) multiple times
# vlengths can either be an int or a list of ints
def montecarlo_sample_tpds(distmatrix, vlengths, graph, prop="goterms", m=1000000,
                           ncores=1, outfile=None, sampling="weighted", approx_after=50000,
                           plot_distributions=False, plot_prefix="distribution_comparison", plot_dir="KDEplots/"):
    start = timeit.default_timer()

    if isinstance(vlengths, int):
        vlengths = [vlengths]

    tpddict = {}
    for i in vlengths:
        tpddict[i] = []

    # get weighted list(vector) of nodes for sampling from
    if sampling == "weighted":  # do weighted sampling
        svec = get_sampling_vector(graph, prop=prop)
    else:  # don't do weighted sampling
        # load graph
        if isinstance(graph, str):
            g = load_graph(graph)
        elif isinstance(graph, Graph):
            g = graph
        else:
            print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        svec = list(g.get_vertices())

    # load dist array from file
    distarray = load_dist_matrix(distmatrix)

    print("array created")
    stop = timeit.default_timer()
    print('Time: ', stop - start)

    for numv in vlengths:
        pbardescr = str("k = " + str(numv))
        pbar = tqdm(total=m, desc=pbardescr)
        actual_samples = min(approx_after, m)
        with concurrent.futures.ThreadPoolExecutor(max_workers=ncores) as executor:
            futuretpd = {executor.submit(get_random_tpd, distarray, svec, numv): i for i in range(actual_samples)}
            for future in concurrent.futures.as_completed(futuretpd):
                tpddict[numv].append(future.result())
                pbar.update()

        # Coerce to numeric and clean data
        coerced_clean_data = []
        for val in tpddict[numv]:
            try:
                numeric_val = float(val)
                if np.isfinite(numeric_val):
                    coerced_clean_data.append(numeric_val)
            except (ValueError, TypeError):
                continue

        if len(coerced_clean_data) < len(tpddict[numv]):
            print(f"Warning: Removed {len(tpddict[numv]) - len(coerced_clean_data)} invalid values from the data.")

        # Check if enough data remains for KDE
        if len(coerced_clean_data) < 2:
            print(f"Not enough data for KDE. Skipping approximation for {numv} vertices.")
            continue

        # Approximation phase
        approx_samples = []
        if m > approx_after:
            # split approximated data into 80:20 train/test set to ensure robust accuracy of KS and plots
            trainsize = int(len(coerced_clean_data)*0.8)
            train_kde_data = coerced_clean_data[0:trainsize]
            test_kde_data = coerced_clean_data[trainsize:]
            #print(train_kde_data)
            #print(type(train_kde_data[1]))
            kde = gaussian_kde(train_kde_data)
            approx_samples = kde.resample(m - approx_after).flatten().tolist()
            tpddict[numv].extend(approx_samples)
            dat_range = abs(max(coerced_clean_data) - min(coerced_clean_data))  # range of data for use in plotting
            dat_buffer = dat_range * 0.05  # a buffer to be added to either side of distribution in pdf plotting
            dat_space = np.linspace(min(coerced_clean_data) - dat_buffer, max(coerced_clean_data) + dat_buffer, num=100)  # total space of plot

            kde_pdf = kde.evaluate(dat_space)  # get pdf estimate of KDE

            # Plotting phase
            if plot_distributions:
                plt.figure(figsize=(10, 6))
                plt.hist(tpddict[numv][:actual_samples], bins=50, alpha=0.5, density=True, label='Sampled Distribution (test split from 80:20 train/test)')
                plt.plot(dat_space, kde_pdf, label='KDE PDF')
                plt.title(f'Distribution Comparison for {numv} Vertices')
                plt.xlabel('Total Pairwise Distance')
                plt.ylabel('Probability')
                plt.legend()

                # Compute KS statistic
                ks_stat, ks_pval = ks_2samp(test_kde_data, kde_pdf)  # KS test of test data versus generated KDE
                plt.figtext(0.5, 0.01, f"KS Statistic: {ks_stat:.4f}, P-Value: {ks_pval:.4f}", ha="center", fontsize=12)


            # plotting comparison of fit between KDE and sampled distributions
            if not plot_dir.endswith(r"/"):  # appends / if necessary
                plot_dir.append(r"/")

            if not os.path.exists(plot_dir):  # checks to see if dir exists, and creates it if it does not
                os.makedirs(plot_dir)
                print(f"Directory '{plot_dir}' was created.")
            else:
                print(f"Directory '{plot_dir}' already exists. Saving KDE plots to this directory.")

            plot_outfile = plot_dir + plot_prefix + "_" + str(numv) + "nodeKDEvsSampledComparison.png"
            plt.savefig(plot_outfile)
            plt.close()

        # Output results
        printstring = ','.join(map(str, tpddict[numv]))
        if outfile is None:
            print("%s,%s\n" % (numv, printstring))
        else:
            with open(outfile, "a+") as mcf:
                mcf.write("%s,%s\n" % (numv, printstring))

        stop = timeit.default_timer()
        #print(str(numv) + " completed after " + str(stop - start) + " seconds")

    stop = timeit.default_timer()
    print('End time: ', stop - start)
    #print(tpddict)

    return tpddict


# loads distance matrix from file into memory, for use by get_tpd, get_random_tpd, etc.
def load_dist_matrix(distmatrix, indexcol=False):
    start = timeit.default_timer()  # for timing duh

    # open dist matrix
    with open(distmatrix, "r") as f:
        if indexcol:  # old cpp script generated an index column but the new method doesn't, so use this arg if needed
            arraysize = len(f.readline().split()) - 1  # array size is -1 because of index column
        else:
            arraysize = len(f.readline().split())
        distarray = np.zeros((arraysize, arraysize))
        flist = f.readlines()
        del flist[0]
        print("file read")
        stop = timeit.default_timer()  # Timer Block
        print('Time: ', stop - start)  # Timer Block

        # store in memory as an array for speedy access
        for line in range(len(flist)):
            try:
                distarray[line] = np.asarray(flist[0].split())[1:]  # add split line to ndarray, skipping first row/col
                del flist[0]  # clear up memory
            except IndexError as ie:
                print(ie)
                pass
        print("array created")
        stop = timeit.default_timer()  # Timer Block
        print('Time: ', stop - start)  # Timer Block

    return distarray


# utility function to calculcate All Pairs Shortest Paths via graph-tool and save all shortest distances as annotations
def annotate_apsp(graphfile, outfile="distAnnotatedGraph.graphml"):
    g = load_graph(graphfile)  # takes first argument as input graph filename
    dist = shortest_distance(g, weights=g.ep.weight)  # calculate all shortest paths and store as vertex property map
    g.vp["shortest_distances"] = dist  # make property internal
    g.save(outfile)


# creates a distance matrix from the graph-annotated distances calculated with gt.shortest_distance()
# row and column indices are based on vertex indices from the graph
def distmatrix_from_annotations(graphfile, outfile="distmatrix.tsv"):
    g = load_graph(graphfile)  # load graph
    arraylist = []  # initialize empty list for 1d arrays to be appended to

    for v in g.vertices():  # iterate through all vertices
        try:
            distvec = np.array(g.vp.shortest_distances[v])
            arraylist.append(distvec)
        except:
            print("no annotation called 'shortest_distances' found. Please use gt.shortest_distances to find APSP")

    distarray = np.vstack(arraylist)  # create 2D array from the list of 1D arrays
    pd.DataFrame(distarray).to_csv(outfile, sep="\t")
    return(distarray)


# get tpd for a single set
def get_tpd(nodelist, distarray):

    combs = combinations(nodelist, 2)  # get all possible combinations for nodes annotated with term

    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory
    row = list(row)
    col = list(col)
    tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
    #print("TPD of set: " + str(tpd))

    return tpd


# get tptpd for a single set
def get_tptpd(nodelist, distmatrix, tpthresh=0.4):
    start = timeit.default_timer()  # Timer Block
    # open dist matrix
    with open(distmatrix, "r") as f:
        arraysize = len(f.readline().split()) - 1  # array size is -1 because of index column
        distarray = np.zeros((arraysize, arraysize))
        flist = f.readlines()
        del flist[0]
        print("file read")
        stop = timeit.default_timer()  # Timer Block
        print('Time: ', stop - start)  # Timer Block

        # store in memory as an array for speedy access
        for line in range(len(flist)):
            try:
                distarray[line] = np.asarray(flist[0].split())[1:]  # add split line to ndarray, skipping first row/col
                del flist[0]  # clear up memory
            except IndexError as ie:
                print(ie)
                pass
        print("array created")

        # do TPTPD calculation
        combs = combinations(nodelist, 2)  # get all possible combinations for nodes annotated with term

        col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
        del combs  # delete to clear up memory

        # get array of combinations with respective shortest path distances
        distlist = distarray[row, col]
        tosumarray = np.array([[[row, col]],
                               [[distlist,
                                 distlist]]])  # create a 2x2xncombs np array with axes source x target x dist
        del distlist

        # for each node, get it's individual tpd and add it to a list to be sorted and trimmed
        tpdlist = []
        for i in range(len(nodelist)):
            mask = tosumarray[0, :, :] == i
            tpd = np.sum(tosumarray[1, :, :],
                         where=mask)  # vectorized addition of PD at all combs for TPD calculation
            tpdlist.append(tpd)

        # add nodes with tpds into data frame, sort it, calculate regular tpd for top% of nodes
        tosortframe = pd.DataFrame({"node": nodelist, "tpd": tpdlist})
        sortedframe = tosortframe.sort_values('tpd', ascending=False)

        numtop = int(len(nodelist) * tpthresh)
        if numtop <= 3:  # make sure there are always at least 3 nodes in set
            numtop = 3
            # TODO make this so it can handle 2?

        newnodeslist = sortedframe['node'].tolist()
        newnodeslist = newnodeslist[0:numtop]  # cut down list to only top nodes

        #  recalculate tpd for only top nodes
        combs = combinations(newnodeslist, 2)  # get all possible combinations for sampled nodes
        col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
        del combs  # delete to clear up memory
        row = list(row)
        col = list(col)
        tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
        #print("Top% (" + str(tpthresh * 100) + "%) TPD: " + str(tpd))

    return tpd


# calculate tpd for all goterms in a given godict (created by make_go_dict())
def all_go_tpd(godict, distmatrix, outfile):
    start = timeit.default_timer()  # for timing duh
    keys = list(godict.keys())  # initialize list of keys

    distarray = load_dist_matrix(distmatrix)
    print("Distance matrix loaded...")

    # calculate tpd for all go terms as found in keys
    for term in tqdm(keys):
        combs = combinations(godict[term], 2)  # get all possible combinations for nodes annotated with term

        try:
            col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
            del combs  # delete to clear up memory
            row = list(row)
            col = list(col)
            tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
            #print(term + "TPD: " + str(tpd))
            with open(outfile, "a+") as f:
                f.write(term + ", " + str(tpd) + "\n")
                #print(term + " TPD written to " + outfile)
        except ValueError as ve:
            print("ERROR ON " + str(term))
            print(ve)
            print("continuing...")
            continue

    stop = timeit.default_timer()  # Timer Block
    print('Total Runtime: ', stop - start)  # Timer Block

    return None


# calculate tpd for all goterms in a given godict (created by make_go_dict())
def all_go_tptpd(godict, distmatrix, outfile, tpthresh=0.2):
    start = timeit.default_timer()  # for timing duh
    keys = list(godict.keys())  # initialize list of keys

    # open dist matrix
    distarray = load_dist_matrix(distmatrix)
    print("Distance matrix loaded...")

    # calculate tpd for all go terms as found in keys
    for term in tqdm(keys):
        combs = combinations(godict[term], 2)  # get all possible combinations for nodes annotated with term

        try:
            col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
            del combs  # delete to clear up memory

            # get array of combinations with respective shortest path distances
            distlist = distarray[row, col]
            tosumarray = np.array([[[row, col]],
                                   [[distlist,
                                     distlist]]])  # create a 2x2xncombs np array with axes source x target x dist
            del distlist

            # for each node, get it's individual tpd and add it to a list to be sorted and trimmed
            tpdlist = []
            for i in range(len(godict[term])):
                mask = tosumarray[0, :, :] == i
                tpd = np.sum(tosumarray[1, :, :],
                             where=mask)  # vectorized addition of PD at all combs for TPD calculation
                tpdlist.append(tpd)

            # add nodes with tpds into data frame, sort it, calculate regular tpd for top% of nodes
            tosortframe = pd.DataFrame({"node": godict[term], "tpd": tpdlist})
            sortedframe = tosortframe.sort_values('tpd', ascending=False)

            numtop = int(len(godict[term]) * tpthresh)
            if numtop <= 3:  # make sure there are always at least 3 nodes in set
                numtop = 3


            newnodeslist = sortedframe['node'].tolist()
            newnodeslist = newnodeslist[0:numtop]  # cut down list to only top nodes

            #  recalculate tpd for only top nodes
            combs = combinations(newnodeslist, 2)  # get all possible combinations for sampled nodes
            col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
            del combs  # delete to clear up memory
            row = list(row)
            col = list(col)
            tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
            #print(term + "Top% (" + str(tpthresh * 100) + "%) TPD: " + str(tpd))

            with open(outfile, "a+") as f:
                f.write(term + ", " + str(tpd) + "\n")

        except ValueError as ve:
            print("ERROR ON " + str(term))
            print(ve)
            print("continuing...")
            continue

    stop = timeit.default_timer()  # Timer Block
    print('Total Runtime: ', stop - start)  # Timer Block

    return None


# gets linear distance between nodes/genes from a given set and plot them against their go rank
def get_linear_distances_tpd(resultsfile, gafname, gencodefile, outfile=None, nsamples=1000000, verbose=True, plot=True):
    # create "gene":"chr:startpos-endpos" dictionary
    if verbose:
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
    if verbose:
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
    if verbose:
        print(" getting ranked list of GO terms from results file...")
    rankedgolist = []
    with open(resultsfile, "r") as resf:  # open results file
        for line in resf:  # iterate through file to extract ranked list of go terms
            if not line.startswith(","):  # skip header line of file
                rankedgolist.append(line.split(",")[1].strip())  # append go term to list

    # in order, iterate through ranked list of GO terms
    if verbose:
        print("getting linear distances for all pairs of genes per ranked GO terms...")
    golindistdict = {}
    for term in rankedgolist:  # iterate through ranked list of goterms
        # get all genes associated with that term
        try:
            genelist = gogenesdict[term]
        except KeyError as ke:
            if verbose:
                print("Could not get gene info for " + term)
            golindistdict[term] = -1
            continue
        # get position/loc info for each term in list
        locdict = {}  # dictionary will contain chromosomes as key, list of positions within that chromosome as value
        for g in genelist:
            try:
                locinfo = geneposdict[g]
            except KeyError as ke:
                if verbose:
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
            if len(locdict[
                       chrom]) >= 2:  # only can count pairwise distance if there are at least two genes on chromosome
                uniquelist = list(
                    set(locdict[chrom]))  # ensure all gene entries in list are unique by removing duplicates
                combs = list(combinations(uniquelist, 2))  # get all pairwise combinations
                paircount = paircount + len(combs)  # add pair count for chromosome to total pair count
                for c in combs:  # for each pair, find distance between them
                    g1, g2 = c  # unpack tuple
                    start1 = g1.split("-")[0]  # get start pos for gene 1
                    stop1 = g1.split("-")[1]  # get stop pos for gene 1
                    start2 = g2.split("-")[0]  # get start pos for gene 2
                    stop2 = g2.split("-")[1]  # get stop pos for gene 2
                    if verbose:
                        print(
                            "start1: " + str(int(start1)) + ", stop1: " + str(int(stop1)) + ", start2: "
                            + str(int(start2)) + ", stop2: " + str(int(stop2)))
                    if int(start1) > int(stop2):  # if gene 1 is later pos that gene 2, dist = start1 - stop2
                        pairdist = int(start1) - int(stop2)
                        #print(str(start1) + "-" + str(stop2))
                    elif int(start2) > int(stop1):  # else is opposite
                        pairdist = int(start2) - int(stop1)
                        #print(str(start2) + "-" + str(stop1))
                    else:  # only remaining case should be if start and stop are equal (which should be exceedingly rare)
                        pairdist = 0
                    if verbose:
                        print("pairdist: " + str(pairdist))
                    dist = dist + pairdist  # add distance of pair to total summed distance
                if verbose:
                    print(term + ": " + str(dist) + "bp")

        # get average linear distance for gene set
        try:
            avglindist = float(dist) / float(paircount)
            golindistdict[term] = avglindist
            if verbose:
                print("average linear distance for set: " + str(avglindist))
        except ZeroDivisionError as ze:
            if verbose:
                print("Could not calculate linear distance for term: " + term)
            continue

    # check dict
    if verbose:
        for k in list(golindistdict.keys()):
            print(k + ": " + str(golindistdict[k]))

    # get average linear distance for random pairs of genes (n samples = 1000000)
    if verbose:
        print("sampling random gene pairs to get average linear distance (N=" + str(nsamples) + ")...")
    samplecounter = 0
    distlist = []
    with tqdm(total=nsamples) as pbar:  # make pbar for sampling loop
        while samplecounter < nsamples:
            genepair = random.sample(list(geneposdict.keys()), 2)
            if geneposdict[genepair[0]].split(":")[0] == geneposdict[genepair[1]].split(":")[
                0]:  # if random genes share the same chromosome
                g1, g2 = genepair  # unpack tuple
                # print(genepair)
                # print(geneposdict[genepair[0]])
                # print(geneposdict[genepair[1]])
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

    # get mean and standard deviation for random gene pairs
    mean = sum(distlist) / len(distlist)
    variance = sum([((x - mean) ** 2) for x in distlist]) / len(distlist)
    stdev = variance ** 0.5

    # create a plot of linear distance vs GO term ranking with line for avg linear distance and std dev
    # create lists to plot
    ranking = range(len(rankedgolist))  # data for x-axis
    rankeddistlist = []  # initialize list for y-axis
    for term in rankedgolist:  # append distance for each ranked GO term in order from most to least significant
        try:
            rankeddistlist.append(golindistdict[term])
        except KeyError as ke:
            if verbose:
                print("no linear distance found for " + term)
            rankeddistlist.append(-1)

    # replace -1 distances with None
    ilist = []
    for i in range(len(rankeddistlist)):
        if rankeddistlist[i] != -1:
            ilist.append(i)  # list of indices with actual linear distance values
    reallindistlist = [rankeddistlist[x] for x in ilist]

    if plot:
        print("plotting...")
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
        plt.axhline(y=mean + stdev, color='red', linestyle='--', label=upperlabel)  # stddev upper bound
        plt.axhline(y=mean - stdev, color='red', linestyle='--')  # stddev lower bound

        # Show a legend
        plt.legend()

        # Display the plot
        plt.savefig(outfile)
        if verbose:
            print("plot saved as " + outfile)

        corr, _ = pearsonr(ranking, rankeddistlist)

        if verbose:
            print("Pearson Correlation = " + str(corr))

    if outfile:
        with open(outfile, "w") as f:
            for i in rankeddistlist:
                f.write(str(i) + "\n")

    if verbose:
        print(rankeddistlist)

    return rankeddistlist


# compare average linear distances of nodes in same MCL/Spectral cluster versus average linear distance globally
def do_linear_analysis_mcl(graphfile, resultsfile, gencodefile, outfile=None, nsamples=10000, verbose=True, plot=True,
                           plotfile="plot.csv"):
    # create "gene":"chr:startpos-endpos" dictionary
    if verbose:
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
    if verbose:
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
                if verbose:
                    print(g1 + " or " + g2 + " not in dict. Skipping...")  # TODO remove?
                continue

        # add avg dist to list
        if distlist:
            avgdist = sum(distlist)/len(distlist)  # get average (mean) of distances from all nodes that share chromosomes
        else:
            avgdist = -1  # if distlist is empty, no nodes are on same chrom so set avgdist to -1
        avg_distlist.append(avgdist)

    # get average linear distance for random pairs of genes (n samples = 10000)
    if verbose:
        print("sampling random gene pairs to get average linear distance (N=" + str(nsamples) + ")...")
    samplecounter = 0
    distlist = []

    with tqdm(total=nsamples) as pbar:  # make pbar for sampling loop
        while samplecounter < nsamples:
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

    if plot:
        # plotting

        # create a plot of linear distance vs random ranking with line for avg linear distance and std dev
        # create lists to plot
        if verbose:
            print("plotting...")
        ranking = range(len(avg_distlist))  # data for x-axis

        # replace -1 distances with None
        ilist = []
        for i in range(len(avg_distlist)):
            if avg_distlist[i] != -1:
                ilist.append(i)  # list of indices with actual linear distance values
        reallindistlist = [avg_distlist[x] for x in ilist]

        # get mean and standard deviation for random gene pairs
        mean = sum(reallindistlist) / len(reallindistlist)
        variance = sum([((x - mean) ** 2) for x in reallindistlist]) / len(reallindistlist)
        stdev = variance ** 0.5

        # Create a scatter plot
        plt.figure().set_figwidth(15)
        plt.scatter(ilist, reallindistlist, s=0.5, label='Linear Distance vs significant gene cluster', color='black')

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
        if verbose:
            print("plot saved as " + plotfile)

    if outfile:
        # write avg distances to file
        with open(outfile, "w") as outf:
            for dist in avg_distlist:
                outf.write(str(dist) + "\n")


    return avg_distlist


# get list of ks for a graph, where k is the number of vertices annotated with a particular GO term
def get_vlengths(graph, prop="goterms"):
    godict = make_go_dict(graph, prop=prop)
    vlengthlist = []
    for key in godict.keys():
        if len(godict[key]) < 2:  # we only care about clusters of size two or more
            continue
        vlengthlist.append(len(godict[key]))
    vlengthlist = sorted(list(set(vlengthlist)))
    return vlengthlist


# pairwise swapping of goterms between random nodes for generation of False Discovery Rate
# 1000000000x(for random pair of nodes, swap one value of 'goterms')
# vars refer to go, as in "golist", but could be reactome or any other annotation set as well
# DEPRECATED in favour of swap_all_goterms()
def label_swap(graph, label="goterms", outfile="swappedGraph.gt", m=3000000):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    pbar = tqdm(total=m)  # initialize loading bar

    numnodes = len(g.get_vertices())  # get number of nodes (vertices) in graph
    numswaps = 1  # track number of actual swaps that are made

    while numswaps <= m:  # loop until desired number of swaps is achieved
        nodepair = random.sample(range(numnodes), 2)  # get a random pair of nodes from the graph

        # extract list of goterms/pathways for random node 1
        v1 = g.vertex(nodepair[0])
        if label == "goterms":
            golist1 = list(g.vp.goterms[v1])
        elif label == "reactomepathways":
            golist1 = list(g.vp.reactomepathways[v1])
        elif label == "chromosomes":
            golist1 = list(g.vp.vname[v1])

        # extract list of goterms/pathways for random node 2
        v2 = g.vertex(nodepair[1])
        if label == "goterms":
            golist2 = list(g.vp.goterms[v2])
        elif label == "reactomepathways":
            golist2 = list(g.vp.reactomepathways[v2])
        elif label == "chromosomes":
            golist2 = list(g.vp.vname[v2])

        # check to make sure neither list is empty, if that's true then swap a term in one list for one in the other
        if len(golist1) != 0 and len(golist2) != 0:
            t1 = golist1[random.sample(range(len(golist1)), 1)[0]]
            t2 = golist2[random.sample(range(len(golist2)), 1)[0]]
            ind1 = golist1.index(t1)
            ind2 = golist2.index(t2)
            golist1[ind1] = t2
            golist2[ind2] = t1
            if label == "goterms":
                g.vp.goterms[v1] = golist1
            elif label == "reactomepathways":
                g.vp.reactomepathways[v1] = golist1
            if label == "goterms":
                g.vp.goterms[v2] = golist2
            elif label == "reactomepathways":
                g.vp.reactomepathways[v2] = golist2

            numswaps = numswaps + 1  # count that swap was made
            pbar.update()
        # if either node has no goterms, pass
        else:
            pass

    g.save(outfile, fmt="gt")
    return None


# new implementation of shuffling that guarantees GO terms are annotated to the same number of nodes as originally
def swap_all_goterms(graph, outfile="swappedGraph.gt"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    # build list of terms to sample from and dictionary of node:#annotations
    annot_per_node_dict = {}
    annot_bag = []  # list to be used as a bag randomizer
    print("Filling bag randomizer...")
    for v in tqdm(g.iter_vertices()):
        annot_list = list(g.vp.goterms[v])  # retrieve list of annotations
        annot_per_node_dict[v] = len(annot_list)  # get count of annotations for later
        annot_bag.extend(annot_list)  # add list of annotations to bag

    # get counts of go terms and list of vertices to recieve new annotations
    annot_set = list(set(annot_bag))
    vlist = list(annot_per_node_dict.keys())
    vlistcopy = [x for x in vlist]
    annot_countdict = {}
    for annot in annot_set:  # populate a dictionary where key:value is [annotation]:[number of vertices with annotation]
        annot_countdict[annot] = annot_bag.count(annot)

    new_annot_dict = {}
    # for each annotation, randomly assign to setsize (annot_countdict[annot]) number of vertices
    for annot in annot_set:
        try:
            x = annot_countdict[annot]
            for v in random.sample(vlist, x):
                if v in new_annot_dict.keys():
                    new_annot_dict[v].append(annot)
                else:
                    new_annot_dict[v] = [annot]
                # if vertex has enough annotations, remove it from vlist
                if len(new_annot_dict[v]) >= annot_per_node_dict[v]:
                    vlist.remove(v)
        except ValueError as ve:
            print("Not enough vertices left to sample. Allowing sampling from all previously annotated vertices...")
            for v in random.sample(vlistcopy, x):
                if v in new_annot_dict.keys():
                    new_annot_dict[v].append(annot)
                else:
                    new_annot_dict[v] = [annot]

    # convert dict to vertex property array
    new_goterm_prop = g.new_vertex_property("vector<string>")
    for v in new_annot_dict.keys():
        new_goterm_prop[v] = new_annot_dict[v]

    # make property map internal and save graph
    g.vp.goterms = new_goterm_prop
    g.save(outfile)


# remove nodes with degree greater than degreecutoff from graph
def trim_nodes_by_degree(graph, outfile="trimmed.gt", degreecutoff=5000):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    del_list = []
    for v in g.vertices():
        if v.out_degree() >= 5000:
            del_list.append(v)
    for v in reversed(sorted(del_list)):
        g.remove_vertex(v)
    g.save(outfile, fmt="graphml")


# generates csv with every GO term and the number of nodes annotated with that term
def make_numnodes_dict(graph, prop="goterms"):
    if prop == "goterms":
        mygodict = make_go_dict(graph)
    elif prop == "reactomepathways":
        mygodict = make_go_dict(graph, prop="reactomepathways")
    nodesdict = {}
    for key in mygodict.keys():
        nodesdict[key] = len(mygodict[key])
    return nodesdict


# get pvals for every go term
def get_go_tpd_pvals(graph, tpdfile, shuftpdfile, distrnfile, approxdist=False, prop="goterms"):
    # read in all files
    if prop == "goterms":
        nndict = make_numnodes_dict(graph)  # nndict is goterm:numnodeswithgoterm
    elif prop == "reactomepathways":
        nndict = make_numnodes_dict(graph, prop="reactomepathways")

    # open Goterm,tpd .csv file as goterm:tpd dict
    tpddict = {}
    with open(tpdfile, "r") as f:
        for line in f:
            splitline = line.split(", ")
            tpddict[splitline[0]] = splitline[1][:-2]

    # open goterm:shuffledtpd csv file as goterm:shuffledtpd dict
    shuftpddict = {}
    with open(shuftpdfile, "r") as f:
        for line in f:
            splitline = line.split(", ")
            shuftpddict[splitline[0]] = splitline[1][:-2]

    # opens knodes:MCdistribution csv as knodes:[MCdistribution] dict
    distdict = {}
    with open(distrnfile, "r") as f:
        for line in f:
            splitline = line.split(",")
            distdict[splitline[0]] = sorted(splitline[2:-3])
            distdict[splitline[0]] = [float(i) for i in distdict[splitline[0]]]  # convert to float

    pvalsdict = {}  # initialize goterm:pval dict
    shufpvalsdict = {}  # initialize goterm:shuffledpval dict

    # for each key (goterm) in our files, look up k, then find that MCdistribution
    # for i in tqdm(tpddict.keys()):
    for i in tpddict.keys():
        if str(nndict[i]) not in distdict:  # skip values of k with no MCdistribution
            print(str(nndict[i]) + " not in MC distribution")
            continue
        if i not in shuftpddict.keys():  # make sure term is in both dicts
            print(str(i) + " not found in both TPD dictionaries")
            continue

        # TODO fix this / make sure it works
        if approxdist:  # approximate the tpd distribution function instead of using the empirical one
            shape, location, scale = sps.gamma.fit(distdict[str(nndict[i])])  # unpack tuple of shape and scale params
            pval = sps.gamma.pdf(np.float64(tpddict[i]), shape, location, scale)  # pulls pval for i from fit distrn
            pvalsdict[i] = pval
            # print("TPD: " + str(tpddict[i]))
            # print("realpval: " + str(pvalsdict[i]))
            # repeat for shuffled set
            pval = sps.gamma.pdf(np.float64(shuftpddict[i]), shape, location, scale)  # pulls pval for i from fit distrn
            shufpvalsdict[i] = pval
            # print("shufTPD: " + str(tpddict[i]))
            # print("shufpval: " + str(pvalsdict[i]))

        else:
            # for every value in distribution, count if it is greater than current goterm's TPD
            # this means a low p-value is associated with our TPD being lower than most vals in MCdist
            counter = 0
            for j in distdict[str(nndict[i])]:
                if float(tpddict[i]) >= float(j):
                    counter += 1
            # print("\n")
            # print(i)
            # print(nndict[i])
            # print("TPD: " + str(tpddict[i]))
            pval = counter / len(distdict[str(nndict[i])])
            pvalsdict[i] = pval
            # print("realpval: " + str(pvalsdict[i]))

            # do the same for the shuffled set
            counter = 0
            for j in distdict[str(nndict[i])]:
                if float(shuftpddict[i]) >= float(j):  # TODO make sure is > not <
                    counter += 1
            # print("shufTPD: " + str(shuftpddict[i]))
            shufpval = counter / len(distdict[str(nndict[i])])
            shufpvalsdict[i] = shufpval
            # print("shufpval: " + str(shufpvalsdict[i]))

        #print(str(i) + " pval: " + str(pval) + " / shufpval: " + str(shufpval))
    # coerce to dataframe
    df = pd.DataFrame([nndict, tpddict, pvalsdict, shuftpddict, shufpvalsdict],
                      index=["nnodes", "tpd", "pval", "shuftpd", "shufpval"])
    df = df.dropna(axis=1)

    return df


# get pvals for generic sets from file
# tpdfiles will need to have columns setname, nnodes, clusterscore(tpd)
def get_tpd_pvals(tpdfile, shuftpdfile, distrnfile, approxdist=True):
    # open normal and shuffled clusterscore(tpd) files as setname, nnodes, clusterscore(tpd) df
    df = pd.read_csv(tpdfile, names=["setname", "nnodes", "tpd"], header=0)
    shufdf = pd.read_csv(shuftpdfile, names=["setname", "nnodes", "shuftpd"], header=0)

    # explicitly type all columns so stupid merging will work
    df.setname.astype(str)
    df.nnodes.astype(int)
    df.tpd.astype(float)
    shufdf.setname.astype(str)
    shufdf.nnodes.astype(int)
    shufdf.shuftpd.astype(float)

    df = df.merge(shufdf, on="setname", suffixes=["DROP", None]).filter(regex="^(?!.*DROP)")  # left joins both dfs and drops repeat cols

    # opens knodes:MCdistribution csv as knodes:[MCdistribution] dict
    distdict = {}
    with open(distrnfile, "r") as f:
        for line in f:
            splitline = line.split(",")
            distdict[splitline[0]] = sorted(splitline[2:-3])
            distdict[splitline[0]] = [float(i) for i in distdict[splitline[0]]]  # convert to float

    pvalsdict = {}  # initialize set:pval dict
    shufpvalsdict = {}  # initialize set.:shuffledpval dict

    # port legend
    # i = row[0]
    # tpddict[i] = row[1]
    # nndict[i] = row[2]
    # shuftpddict[i] = row[3]

    # for each set in our files find that MCdistribution and get pval
    for row in df.itertuples(index=False):  # access values like row[0]
        if str(row[2]) not in distdict:  # skip values of k with no MCdistribution
            continue

        if approxdist:  # approximate the tpd distribution function instead of using the empirical one
            shape, location, scale = sps.gamma.fit(
                distdict[str(row[2])])  # unpack tuple of shape and scale params
            pval = sps.gamma.pdf(np.float64(row[1]), shape, location, scale)  # pulls pval for i from fit distrn
            pvalsdict[row[0]] = pval
            # print("TPD: " + str(tpddict[i]))
            # print("realpval: " + str(pvalsdict[i]))
            # repeat for shuffled set
            pval = sps.gamma.pdf(np.float64(row[3]), shape, location, scale)  # pulls pval for i from fit distrn
            shufpvalsdict[row[0]] = pval
            # print("shufTPD: " + str(tpddict[i]))
            # print("shufpval: " + str(pvalsdict[i]))

        else:
            # for every value in distribution, count if it is greater than current goterm's TPD
            # this means a low p-value is associated with our TPD being lower than most vals in MCdist
            counter = 0
            for j in distdict[str(row[2])]:
                if float(row[1]) >= float(j):
                    counter += 1
            # print("\n")
            # print(i)
            # print(nndict[i])
            # print("TPD: " + str(tpddict[i]))
            pval = counter / len(distdict[str(row[2])])
            pvalsdict[row[0]] = pval
            # print("realpval: " + str(pvalsdict[i]))

            # do the same for the shuffled set
            counter = 0
            for j in distdict[str(row[2])]:
                if float(row[3]) >= float(j):
                    counter += 1
            # print("shufTPD: " + str(shuftpddict[i]))
            shufpval = counter / len(distdict[str(row[2])])
            shufpvalsdict[row[0]] = shufpval
            # print("shufpval: " + str(shufpvalsdict[i]))

    # coerce to dataframe
    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
    pvalsdf = pd.DataFrame.from_dict(pvalsdict, orient="index", columns=["pval"])
    pvalsdf.reset_index(inplace=True)
    pvalsdf = pvalsdf.rename(columns={'index': 'setname'})
    pvalsdf['setname'] = pvalsdf['setname'].astype(str)

    shufpvalsdf = pd.DataFrame.from_dict(shufpvalsdict, orient="index", columns=["shufpval"])
    shufpvalsdf.reset_index(inplace=True)
    shufpvalsdf = shufpvalsdf.rename(columns={'index': 'setname'})
    shufpvalsdf['setname'] = shufpvalsdf['setname'].astype(str)

    df['setname'] = df['setname'].astype(str)

    df = df.merge(pvalsdf, on="setname", suffixes=["DROP", None]).filter(regex="^(?!.*DROP)")
    df = df.merge(shufpvalsdf, on="setname", suffixes=["DROP", None]).filter(regex="^(?!.*DROP)")
    df = df.dropna(axis=1)

    return df


# system call to c++ program that calculates distance matrix
def generate_distance_matrix(graphmlfile, outfile="distmatrix.tsv", sofile="gothicAPSP.so", processors=1):
    # check that files have leading "./" so subprocess.run() will recognize them
    if graphmlfile[0] != ".":
        graphmlfile = "./" + graphmlfile
    if outfile[0] != ".":
        outfile = "./" + outfile
    if sofile[0] != ".":
        sofile = "./" + sofile
    # check if file is graphml or gt; if gt convert to graphml
    if graphmlfile.split(".")[-1] == "gt":
        graphmlname = str(".".join(graphmlfile.split(".")[:-1])) + ".graphml"  # so "."s in name don't cause truncated new name
        g = load_graph(graphmlfile)
        g.save(graphmlname)
        graphmlfile = graphmlname  # replaces argument var with new file
    try:
        cpp_process = subprocess.run([sofile, graphmlfile, outfile], capture_output=True, check=True)
        print("stdout:", cpp_process.stdout)
        print("stderr:", cpp_process.stderr)
    except subprocess.CalledProcessError:
        print("ERROR: Djikstra APSP c++ subprocess failed\n")
        print("stdout:", cpp_process.stdout)
        print("stderr:", cpp_process.stderr)


def count_edges(g):
    count = 0
    for edge in g.edges():
        count = count + 1
    return count


def count_vertices(g):
    count = 0
    for v in g.vertices():
        count = count + 1
    return count


# extracts chromosome from graphmlID for each node and applies it as a new vertex property map
def chromosome_annotate(g):
    vprop_chromosome = g.new_vertex_property("string")

    for v in g.vertices():
        node = int(v)
        vprop_chromosome[node] = g.vp.vname[v].split(":")[0]

    g.vp["chromosome"] = vprop_chromosome  # make vprop map internal to g so it will save with it


# annotate whether edges repressent intrachromosomal or interchromosomal contacts
def contact_annotate(g):
    eprop_contactType = g.new_edge_property("bool")

    try:
        # annotate graph with chromosomes if not already
        if g.vertex_properties["chromosome"]:
            pass

    except KeyError:
        chromosome_annotate(g)

    # test and annotate edges
    for s, t in g.iter_edges():
        if g.vp.chromosome[s] == g.vp.chromosome[t]:
            eprop_contactType[g.edge(s, t)] = False
        else:
            eprop_contactType[g.edge(s, t)] = True

    # make edge property internal
    g.ep["interchromosomal"] = eprop_contactType  # make eprop map internal so it saves with graph


# extract interchromosomal contacts
def extract_interchromosomal_subgraph(graph, outfile, verbose=False):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    # annotate contacts
    contact_annotate(g)

    # save subgraph with only interchromosomal edges
    # create graph view, and then save the graph
    g.set_edge_filter(g.ep.interchromosomal)
    fg = Graph(g=g, directed=False, prune=True)
    gv = extract_largest_component(fg, prune=True)
    if verbose:
        print("New NumEdges: " + str(count_edges(gv)))
    gv.save(outfile)

    if verbose:
        lccCount = count_vertices(gv)
        ogCount = count_vertices(g)
        print("number of vertices in original graph: " + str(ogCount))
        print("number of vertices in largest component of subgraph: " + str(lccCount))

# creates a graph with low-weight edges removed
# if weights_are_distances is True, takes edges with lowest values. If False takes highest values
def create_top_edges_graph(graph, threshold=0.05, outfile="top5percentGraph.gt", weights_are_distances=True):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys. exit()

    print("NumEdges: " + str(count_edges(g)))
    # for all edges in graph, add weight to a list, sort the list, then calculate threshold for top%
    if weights_are_distances:
        eWeightList = np.sort(g.ep.weight.a)  # sorts in ascending order if weights are like distances
    else:
        eWeightList = np.sort(g.ep.weight.a)[::-1]  # sorts in descending order if weights are like similarities

    cutoffIndex = int(len(eWeightList) * threshold)
    cutoffVal = eWeightList[cutoffIndex]
    print("sorted")
    print(cutoffVal)

    # apply a boolean property map where prop=T is edge weight is in the top ?%
    TopPercentEdge = g.new_edge_property("bool")
    g.edge_properties["TopPercentEdge"] = TopPercentEdge

    if weights_are_distances:
        for e in g.edges():
            if g.ep.weight[e] <= cutoffVal:  # weights with distance shorter than cutoff are kept
                g.ep.TopPercentEdge[e] = True
            else:
                g.ep.TopPercentEdge[e] = False
    else:
        for e in g.edges():
            if g.ep.weight[e] >= cutoffVal:  # weights with higher similarity than cutoff are kept
                g.ep.TopPercentEdge[e] = True
            else:
                g.ep.TopPercentEdge[e] = False

    # create graph view, and then save the graph
    g.set_edge_filter(TopPercentEdge)
    fg = Graph(g=g, directed=False, prune=True)
    print("New NumEdges: " + str(count_edges(fg)))
    fg.save(outfile)


# creates a graph with low-weight edges removed
def create_top_edges_plus_mst_graph(graph, threshold=0.05, outfile="top5percentGraph.gt"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys. exit()

    print("NumEdges: " + str(count_edges(g)))
    # for all edges in graph, add weight to a list, sort the list, then calculate threshold for top%
    eWeightList = np.sort(g.ep.weight.a) # sorts in ascending order
    cutoffIndex = int(len(eWeightList) * threshold)  # gets negative index bc list is sorted smallest to largest
    cutoffVal = eWeightList[cutoffIndex]
    print("sorted")

    # get minimum spanning tree
    mst = min_spanning_tree(g, weights=g.ep.weight)  #remember bool(mst[g.edge(s,t)])

    # apply a boolean property map where prop=T is edge weight is in the top ?%
    TopPercentEdge = g.new_edge_property("bool")
    g.edge_properties["TopPercentEdge"] = TopPercentEdge
    for e in g.edges():
        if g.ep.weight[e] <= cutoffVal or bool(mst[e]):
            g.ep.TopPercentEdge[e] = True
        else:
            g.ep.TopPercentEdge[e] = False

    # create graph view, and then save the graph
    g.set_edge_filter(TopPercentEdge)
    print("New NumEdges: " + str(count_edges(g)))
    g.save(outfile)


# takes a graph, removes all but the top k edges per node
def retain_top_k_edges_per_node(graph, k=4, outfile=None):
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
    if not g.is_directed():  # only works if graph is undirected

        # Get the edge weights into a numpy array
        weights = np.array([e for e in g.ep.weight], dtype=float)
        # Get arrays of source and target indices for each edge
        source = np.array([e.source() for e in g.edges()], dtype=int)
        target = np.array([e.target() for e in g.edges()], dtype=int)

        # Combine the source, target, and weights for sorting and manipulation
        edges = np.vstack((source, target, weights)).T

        # For undirected graphs, we treat each edge twice: once for each direction
        edges = np.vstack((edges, edges[:, [1, 0, 2]]))

        # Sort edges based on weights
        edges_sorted = edges[edges[:, 2].argsort()]

        # Identify unique edges after sorting to avoid double removal
        _, unique_indices = np.unique(edges_sorted[:, :2], axis=0, return_index=True)
        edges_unique_sorted = edges_sorted[np.sort(unique_indices)]

        # This step involves custom logic to keep the top k edges for each node
        to_keep_indices = np.array([], dtype=int)
        for node in np.unique(edges_unique_sorted[:, :2]):
            # Find indices of edges connected to this node
            idx = np.where((edges_unique_sorted[:, 0] == node) | (edges_unique_sorted[:, 1] == node))[0]
            # If there are more than k edges, keep only the indices of the top k based on weight
            if len(idx) > k:
                to_keep_indices = np.concatenate((to_keep_indices, idx[:k]))
            else:
                to_keep_indices = np.concatenate((to_keep_indices, idx))

        # Determine edges to remove
        all_indices = np.arange(edges_unique_sorted.shape[0])
        to_remove = np.setdiff1d(all_indices, to_keep_indices)

        # Convert to list of edge descriptors for removal
        edges_to_remove = [g.edge(g.vertex(int(edges_unique_sorted[i, 0])), g.vertex(int(edges_unique_sorted[i, 1])))
                           for i in to_remove]
        # Convert to list of tuples for easier duplicate removal
        edges_to_remove_tuples = [(int(edges_unique_sorted[i, 0]), int(edges_unique_sorted[i, 1])) for i in to_remove]
        # Create a set to hold unique edges, treating each edge as a sorted tuple to handle undirected nature
        unique_edges_to_remove = set(tuple(sorted(edge)) for edge in edges_to_remove_tuples)
        # Now convert these unique, sorted tuples back into edge descriptors for removal
        unique_edges_to_remove_descriptors = [g.edge(g.vertex(edge[0]), g.vertex(edge[1])) for edge in
                                              unique_edges_to_remove]
        # Remove the edges
        for e in unique_edges_to_remove_descriptors:
            if e:  # Check if the edge descriptor is not None
                g.remove_edge(e)

        if outfile is None:
            return g
        else:
            g.save(outfile)

    else:
        raise ValueError("This function is intended for undirected graphs only.")


# takes a graph, removes all but the Minimum Spanning Tree and top k edges per node
def retain_top_k_plus_mst_edges_per_node(graph, k=0, outfile=None):
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
    if not g.is_directed():  # only works if graph is undirected

        # get minimum spanning tree
        mst = min_spanning_tree(g, weights=g.ep.weight)

        # if keeping more than the MST, calculate graph containing top edges per node
        if k > 0:
            # Get the edge weights into a numpy array
            weights = np.array([e for e in g.ep.weight], dtype=float)
            # Get arrays of source and target indices for each edge
            source = np.array([e.source() for e in g.edges()], dtype=int)
            target = np.array([e.target() for e in g.edges()], dtype=int)

            # Combine the source, target, and weights for sorting and manipulation
            edges = np.vstack((source, target, weights)).T

            # For undirected graphs, we treat each edge twice: once for each direction
            edges = np.vstack((edges, edges[:, [1, 0, 2]]))

            # Sort edges based on weights
            edges_sorted = edges[edges[:, 2].argsort()]

            # Identify unique edges after sorting to avoid double removal
            _, unique_indices = np.unique(edges_sorted[:, :2], axis=0, return_index=True)
            edges_unique_sorted = edges_sorted[np.sort(unique_indices)]

            # This step involves custom logic to keep the top k edges for each node
            to_keep_indices = np.array([], dtype=int)
            for node in np.unique(edges_unique_sorted[:, :2]):
                # Find indices of edges connected to this node
                idx = np.where((edges_unique_sorted[:, 0] == node) | (edges_unique_sorted[:, 1] == node))[0]
                # If there are more than k edges, keep only the indices of the top k based on weight
                if len(idx) > k:
                    to_keep_indices = np.concatenate((to_keep_indices, idx[:k]))
                else:
                    to_keep_indices = np.concatenate((to_keep_indices, idx))

            # Determine edges to remove
            all_indices = np.arange(edges_unique_sorted.shape[0])
            to_remove = np.setdiff1d(all_indices, to_keep_indices)

            # Convert to list of edge descriptors for removal
            edges_to_remove = [g.edge(g.vertex(int(edges_unique_sorted[i, 0])), g.vertex(int(edges_unique_sorted[i, 1])))
                               for i in to_remove]
            # Convert to list of tuples for easier duplicate removal
            edges_to_remove_tuples = [(int(edges_unique_sorted[i, 0]), int(edges_unique_sorted[i, 1])) for i in to_remove]
            # Create a set to hold unique edges, treating each edge as a sorted tuple to handle undirected nature
            unique_edges_to_remove = set(tuple(sorted(edge)) for edge in edges_to_remove_tuples)
            # Now convert these unique, sorted tuples back into edge descriptors for removal
            unique_edges_to_remove_descriptors = [g.edge(g.vertex(edge[0]), g.vertex(edge[1])) for edge in
                                                  unique_edges_to_remove]
            # Remove the edges
            for e in unique_edges_to_remove_descriptors:
                if e:  # Check if the edge descriptor is not None
                    if not bool(mst[e]):
                        g.remove_edge(e)

        else:  # if k=0, remove all edges not belonging to the mst
            elist = g.get_edges()
            for s,t in elist:
                if g.edge(s,t):
                    if not bool(mst[g.edge(s,t)]):
                        #print(f"({s},{t}) is {bool(mst[g.edge(s,t)])} removing ({s},{t})")  #TODO remove
                        g.remove_edge(g.edge(s,t))

        if outfile is None:
            return g
        else:
            g.save(outfile)

    else:
        raise ValueError("This function is intended for undirected graphs only.")


# annotates HiC graph with chromosomes, then creates a subgraph for each (only including intra-chrom contacts)
def split_into_chromosomes(graphfile):
    g = load_graph(graphfile)

    try:
        # annotate graph with chromosomes if not already
        if g.edge_properties["contactType"]:
            pass

    except KeyError:
        contact_annotate(g)

    # chromosome list for checking chroms
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    # create a subgraph for each chromosome, and save that subgraph
    for chrom in chrlist:
        # create new temporary edge and vertex properties to store whether v or e is correct chromosome
        vprop_isChromosome = g.new_vertex_property("bool")
        eprop_isChromosome = g.new_edge_property("bool")

        # loops to mark each edge and chromosome
        for v in g.vertices():
            if g.vp.chromosome[v] == chrom:
                vprop_isChromosome[v] = True
            else:
                vprop_isChromosome[v] = False

        for e in g.edges():
            s = e.source()
            if g.vp.chromosome[s] == chrom and g.ep.contactType[e] == "intrachromosomal":
                eprop_isChromosome[e] = True
            else:
                eprop_isChromosome[e] = False

        # set filters
        g.set_filters(eprop_isChromosome, vprop_isChromosome)

        # save graph and unset filters
        graphname = graphfile.split(".")[0] + "_" + str(chrom) + "." + graphfile.split(".")[1]
        g.save(graphname)


# converts bin file (.bed format, could be TADs) to dictionary where key is chr# and value is list of bin bounds
def bins_to_dict(binfile):

    bindict = {"chr1":[], "chr2":[], "chr3":[], "chr4":[], "chr5":[], "chr6":[], "chr7":[], "chr8":[], "chr9":[],
               "chr10":[], "chr11":[], "chr12":[], "chr13":[], "chr14":[], "chr15":[], "chr16":[], "chr17":[],
               "chr18":[], "chr19":[], "chr20":[], "chr21":[], "chr22":[], "chrX":[], "chrY":[], "chrM":[]}

    f = open(binfile, "r")
    for line in f:
        chrom, startpos, endpos = line.split()
        binbounds = str(startpos) + "-" + str(endpos)
        bindict[chrom].append(binbounds)

    return bindict


def generate_graph_report(graph, outfile="graphReport.txt", use_lcc=True):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    chromosome_annotate(g)
    contact_annotate(g)

    # create dict to store all stats to be printed/ used later
    statsdict = {}

    # file size
    if type(graph) == str:
        statsdict["filesize"] = os.path.getsize(graph)

    # number of nodes
    nnodes = count_vertices(g)
    statsdict["nnodes"] = nnodes

    # number of edges
    statsdict["nedges"] = count_edges(g)

    # # if use_lcc is True, use only the largest connected component from now on
    # if use_lcc:
    #     og = g
    #     g = extract_largest_component(og)
    #     print("using largest connected component...")
    #     print(str(count_vertices(g)) + "/" + str(count_vertices(og)) + " nodes in lcc...")

    # degree distribution
    ddistCounts, ddistBins = vertex_hist(g, "total")
    statsdict["ddistCounts"] = ddistCounts
    statsdict["ddistBins"] = ddistBins[1:]  # remove first item so this list repressents the upper range of each bin

    # betweenness distribution (too mem/time intensive for a quick report? needs to calculate APSP. Can use pivots to estimate)
    # extract a list of betweenness values from the generated Vertex Property Map so they can be histogrammed later
    sampleSize = int(nnodes / 10)  # set sample size to be 1/10th of the population
    # get list of unmasked vertices (we have to do this as a loop because graph tool is stupid
    realvlist = []
    try:
        for v in g.vertices():
            realvlist.append(v)
        pivotsList = random.sample(realvlist, sampleSize)  # randomly get a sample of nodes to act as pivots
        vBetweenness, eBetweenness = betweenness(g, pivots=pivotsList, weight=g.ep.weight)
        vBetweenList = []
        for key in vBetweenness:
            vBetweenList.append(float(vBetweenness[key]))
        statsdict["EstBetweennessDist"] = np.array(vBetweenList)  # array so it plays nicely with gnuplot (hopefully)
    except IndexError as ie:
        print(ie)
        print("Cannot get betweenness centrality. Continuing...")
    except ValueError as ve:
        print(ve)
        print("Cannot get betweenness centrality. Continuing...")
    except OverflowError as oe:
        print(oe)
        print("Cannot get betweenness centrality. Continuing...")

    # diameter
    try:
        diameter, dpoints = pseudo_diameter(g, source=realvlist[0], weights=g.ep.weight)
        statsdict["diameter"] = diameter
    except IndexError as ie:
        print(ie)
        print("Cannot get pseudodiameter. Continuing...")

    # global clustering coefficient aka transitivity (weighted)
    try:
        transitivity, transStdDev = global_clustering(g, weight=g.ep.weight)
        statsdict["transitivity"] = transitivity
    except IndexError as ie:
        print(ie)
        print("Cannot get transitivity. Continuing...")

    # number of nodes in largest connected component
    try:
        gv = extract_largest_component(g)  # create a graph view of only the connected component
        lccCount = count_vertices(gv)
        statsdict["nnodesLargestConnectedComponent"] = lccCount
    except ValueError as ve:
        print(ve)
        print("Cannot retrieve largest component. Continuing...")

    # inter vs intra chromosomal contacts ratio
    interCount = 0
    intraCount = 0
    for e in g.edges():
        try:  # catch if graph hasn't been annotated with chromosome / contact info
            if g.ep.interchromosomal[e] == True:
                interCount = interCount + 1
            elif g.ep.interchromosomal[e] == False:
                intraCount = intraCount + 1
        except AttributeError as ae:  # if anything goes wrong, try annotating before counting
            print(ae)
            print("Contact info not found. Annotating...")
            chromosome_annotate(g)
            contact_annotate(g)
            if g.ep.interchromosomal[e] == True:
                interCount = interCount + 1
            elif g.ep.interchromosomal[e] == False:
                intraCount = intraCount + 1

    statsdict["interContactCount"] = interCount
    statsdict["intraContactCount"] = intraCount


    # Generate Report in .txt format for easy readability in terminal
    f = open(outfile, "w")  # open file for writing

    # create list of regular lines to becomethe body of the report
    linelist = []

    # draw ascii logo
    linelist.append(r"           _/         ,          .                           " + "\n")
    linelist.append(r"       , -' )         (\-------.,')            (\________________________________________________," + "\n")
    linelist.append(r"     , ,-/ |          /\) )     \/            ,' _.--------------------------------------------, /" + "\n")
    linelist.append(r"   ,',  /, |         /     >--. ,)           / /\\\        _____ _____ _____ _   _ _____ _____  '" + "\n")
    linelist.append(r"  / ,  //|,'        /'    '\--'\\\)          /,' \\\      |  __ \  _  |_   _| | | |_   _/  __ \ " + "\n")
    linelist.append(r" / ,  // ||       ,'    (.--^( `')         //     \\\     | |  \/ | | | | | | |_| | | | | /  \/" + "\n")
    linelist.append(r"( ,  //  ||,___,-'    (__\\\  '^^^'        //      \\\    | | __| | | | | | |  _  | | | | |    " + "\n")
    linelist.append(r" \  //   ||--.__     (    \`^--)  _____.-'/         \\\   | |_\ \ \_/ / | | | | | |_| |_| \__/\ " + "\n")
    linelist.append(r"  >'/    ||,        (      \|_(\-'      ,'           \\\   \____/\___/  \_/ \_| |_/\___/ \____/ " + "\n")
    linelist.append(r" /,'     ||          \          \      /              \\\                    " + "\n")
    linelist.append(r"(/       ||,         \          )  ,'(     `     `    \\\,    " + "\n\n")

    # draw real data
    # try/except blocks to ensure the report will run despite errors and will be as complete as possible
    now = datetime.datetime.now()
    try:
        if type(graph) == str:
            linelist.append("Report generated for " + graph + " on " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
        else:
            linelist.append("Report generated on " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    except KeyError:
        print("Did not append datetime information because it is missing")
    try:
        linelist.append("Filesize: " + str(statsdict["filesize"] / 1000000000) + "Gb\n")
    except KeyError:
        print("Did not append filesize information because it is missing")
    try:
        linelist.append("Number of Nodes: " + str(statsdict["nnodes"]) + "\n")
    except KeyError:
        print("Did not append node count information because it is missing")
    try:
        linelist.append("Number of Edges: " + str(statsdict["nedges"]) + "\n")
    except KeyError:
        print("Did not append edge count information because it is missing")
    try:
        linelist.append("Number of Nodes in Largest Connected Component: " + str(statsdict["nnodesLargestConnectedComponent"]) + "\n")
    except KeyError:
        print("Did not append largest component information because it is missing")
    try:
        linelist.append("Average node degree: " + str(np.mean(g.get_total_degrees(g.get_vertices()))) + "\n")
    except KeyError:
        print("Did not append node degree information because it is missing")
    try:
        linelist.append("Average edge weight: " + str(np.mean(g.ep.weight.a)) + "\n")
    except KeyError:
        print("Did not append edge weight information because it is missing")
    try:
        linelist.append("Number of Intrachromosomal Contacts: " + str(statsdict["intraContactCount"]) + "\n")
    except KeyError:
        print("Did not append intrachromosomal contact information because it is missing")
    try:
        linelist.append("Number of Interchromosomal Contacts: " + str(statsdict["interContactCount"]) + "\n")
    except KeyError:
        print("Did not append interchromosomal contact information because it is missing")
    try:
        linelist.append("Pseudodiameter: " + str(statsdict["diameter"]) + "\n")
    except KeyError:
        print("Did not append pseudodiameter information because it is missing")
    try:
        linelist.append("Global Clustering Coefficient (Weighted Transitivity): " + str(statsdict["transitivity"]) + "\n")
    except KeyError:
        print("Did not append transitivity information because it is missing")

    # create plots that will be drawn at the bottom of the report
    # plot degree distribution

    max_count = np.max(ddistCounts)  # Determine the maximum count for scaling
    histogram_width = 50  # Width of the histogram in characters
    scale_factor = histogram_width / max_count  # Scale factor to fit counts within the histogram width
    bar_value = 1/scale_factor

    linelist.append("\nDegree Distribution Histogram:\n")
    linelist.append(f"\n# = ~{str(bar_value)}\n\n")

    for i, count in enumerate(ddistCounts):  # Print the histogram
        scaled_count = int(count * scale_factor)  # Scale the count to the current scale factor
        bar = '#' * scaled_count  # Create the bar for the histogram using '#' characters
        linelist.append(f"{ddistBins[i]:>5} | {bar} \n")  # Print the bin and its corresponding bar

    for line in linelist:  # write lines to file
        f.write(line)

    sys.stdout = f
    f.close()
    sys.stdout.close()
    sys.stdout = sys.__stdout__


# function to print list of Top GO terms with FDR and p-val cutoffs for term inclusion
def get_top_goterms(resultscsv, outfile="topterms.csv"):
    # read in and format results csv from getgotpdvals()
    df = pd.read_csv(resultscsv)
    df = df.transpose()
    df.columns = df.iloc[0]
    df["goterm"] = df.index
    df = df[1:]
    df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                    "shufpval": "float", "goterm": "str"})

    sorteddf = df.sort_values("pval")
    shufsorteddf = df.sort_values("shufpval")

    golist = []
    plist = []
    fdrlist = []
    tpdlist = []

    for index, row in sorteddf.iterrows():
        thisgo = row["goterm"]
        thispval = row["pval"]
        thistpd = row["tpd"]

        # get count in list of pval and shufpval <= pval
        pcount = sum(i <= thispval for i in sorteddf["pval"])
        shufcount = sum(i <= thispval for i in shufsorteddf["shufpval"])  # don't even think I needed to sort for this lol
        denom = pcount + shufcount # denominator for fdr calc TODO make sure this is right (not pcount/shufcount?)
        fdr = shufcount / denom

        golist.append(thisgo)
        plist.append(thispval)
        fdrlist.append(fdr)
        tpdlist.append(thistpd)

    # create new df with cols: "goterm, pval, FDR"
    newdf = pd.DataFrame({"goterm": golist, "tpd": tpdlist, "pval": plist, "FDRatThisCutoff": fdrlist})

    # save df
    newdf.to_csv(outfile)


# take results file and return table with cols: pvalThreshold, numRealPassing, numShufPassing, fdr
# applies a monotonic traqnsformation to the FDR such that the FDR never decreases
# ^this is valid because if FDR decreases you can always just pick a higher pval threshold that gave better fdr
# important that startstopstep always decreases for monotonic transformation to work
def get_real_v_null_significants(resultscsv, startstopstep=(0.1, 0, -0.00001),
                                 outfile="SignificantClustersRealVsShuffledTPD.csv", monotonic=True):

    # read in and format results csv from getgotpdvals()
    df = pd.read_csv(resultscsv)
    df = df.transpose()
    df.columns = df.iloc[0]
    df["goterm"] = df.index
    df = df[1:]
    df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                    "shufpval": "float", "goterm": "str"})

    lowestfdr = 1  # initialize lowest fdr as 1, which will be checked against every new fdr value

    pstart, pstop, pstep = startstopstep

    with open(outfile, "w") as f:
        f.write("pvalThreshold, numReal, numShuf, fdr\n")

        for i in np.arange(pstart, pstop, pstep):

            # get count in list of pval and shufpval <= pval
            pcount = sum(j <= i for j in df["pval"])
            shufpcount = sum(j <= i for j in df["shufpval"])
            fdr = shufpcount / pcount

            if monotonic:  # check if monotonic transformation option is selected (True by default)
                if fdr < lowestfdr:  # always ensures fdr that gets written is the lowest possible
                    lowestfdr = fdr

            else:
                lowestfdr = fdr  # if monotonic False, just write new FDR every time

            # write line to file
            f.write(str(i) + ", " + str(pcount) + ", " + str(shufpcount) + ", " + str(lowestfdr) + "\n")


# takes a graph and outputs a list of node ids annotated with a given go term
def get_nodes_with_term(graph, term, cytoscape=False):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    nodelist = []
    for v in g.vertices():
        for t in g.vp.goterms[v]:
            if t == term:
                if cytoscape:
                    nodelist.append("n" + str(int(v)))
                else:
                    nodelist.append(int(v))

    return nodelist


# takes list of clusters and graph, writes one file per cluster with all the genes on that cluster's nodes
def prep_clusters_ontologizer(mclfile, graph, outdir="./mclClusterFiles/", verbose=False, clusterfile_delim=","):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    numlines = sum(1 for line in open(mclfile))  # total number of lines in file
    if verbose:
        print("files to create: " + str(numlines))
    namecounter = 1  # ensures each cluster gets a unique filename equal to the line number

    # check if outdir exists and if not create it
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # load mcl results file
    f = open(mclfile, "r")

    print("working...")

    # for each line in mcl file, take as list of nodes, then iterate through list adding all genes at node to genelist
    for line in f:
        if verbose:
            print(str(namecounter) + "/" + str(numlines), flush=True)
        strlist = line.split(clusterfile_delim)
        nodelist = [int(x) for x in strlist[:-1]]  # convert strings to ints
        if verbose:
            print(nodelist)
        genelist = []
        for node in nodelist:
            if verbose:
                    print(str(node) + ": " + str(g.vp.genes[node]))
            if g.vp.genes[node]:
                genelist.extend(list(g.vp.genes[node]))

            else:
                print(str(node) + " no genes found")

        genelist = list(set(genelist))  # before writing, make sure all genes in list are unique (TODO: right??)
        if verbose:
            print(genelist)

        # open new file in outdir and write list of genes to it
        # if there are no genes in the cluster, append EMPTY to filename
        if not genelist:
            outfilename = outdir + "cluster" + str(namecounter) + "_EMPTY.txt"
        else:
            numgenes = len(genelist)
            outfilename = outdir + "cluster" + str(namecounter) + "_" + str(numgenes) + "genes.txt"

        with open(outfilename, "w") as f2:
            for gene in genelist:
                f2.write(gene + "\n")
            if verbose:
                print(outfilename + " written")
        namecounter = namecounter + 1  # increment name counter

    f.close()


# calculates A/B compartments from .hic file and annotates graph
def annotate_ab_compartments_hic(graphfile, hicfile, outfile, genomefile):
    hic = fanc.load(hicfile)  # load hic file into mem
    ab = fanc.ABCompartmentMatrix.from_hic(hic)  # get
    ev = ab.eigenvector()
    domains = ab.domains()

    g = load_graph(graphfile)

    for region in domains.regions:
        # extract region information
        chr = region.chromosome
        start = region.start
        end = region.end
        compartment = region.name

        print(str(chr) + ":" + str(start) + "-" + str(end) + " -> " + str(compartment))


# annotate nodes with genomic positions using a regions file (.bed)
def annotate_node_genome_positions(graphfile, regionsfile):
    # load graph and regions file to get chr, start and end positions for each node and make that info a graph annotation
    g = load_graph(graphfile)
    with open(regionsfile, "r") as rf:
        pos_vp = g.new_vertex_property("string")
        vcount = 0
        for line in rf:
            chrom, startpos, endpos = line.strip().split()
            pos_vp[vcount] = chrom + ":" + startpos + "-" + endpos
            vcount = vcount + 1
        g.vp.pos = pos_vp


# parse genome fasta according to bins in binbedfile, then count and output G/C ratios for each bin
# the regions file better be sorted in ascending order by chromosome and chrom position or else
def get_node_gc_ratios(binbedfile, genomefastafile, outfile=None):

    bf = open(binbedfile, "r")  #
    gf = open(genomefastafile, "r")  # me and who?

    # iterate through bedfile to get lists of node start positions
    bin_chrom_start_list = []

    for line in bf:
        bin_chrom_start_list.append(int(line.split()[1]))
    bf.close()

    # iterate through genome file and count GC content for each node, at each bin boundary store the GC vs ATGC ratio in a list
    node_gc_ratios_list = []  # list for storing the gc content ratio (NUMGC/NUMALLBASES) of each region
    pos_counter = 0
    base_counter = 0
    gc_counter = 0
    done_flag = 0
    nextpos = bin_chrom_start_list.pop(0)

    for line in gf:

        if done_flag:  # if no more bins to get gc ratios for, break loop
            break

        if line.startswith(">"):
            if base_counter != 0:
                gc_ratio = float(gc_counter)/float(base_counter)
                node_gc_ratios_list.append(gc_ratio)
                pos_counter = 0
                base_counter = 0
                gc_counter = 0

        else:
            base_list = [x for x in line]
            base_list.pop()  # remove newline character from line list so it is not counted as a base
            for b in base_list:
                # add base to counts depending on identity
                base_counter = base_counter + 1
                if b == "G" or b == "C":
                    gc_counter = gc_counter + 1

                pos_counter = pos_counter + 1

                if pos_counter >= nextpos:
                    gc_ratio = float(gc_counter) / float(base_counter)
                    node_gc_ratios_list.append(gc_ratio)
                    #print(str(gc_counter) + "/" + str(base_counter))  #TODO REMOVE WHEN DONE
                    #print(gc_ratio)  #TODO REMOVE WHEN DONE

                    base_counter = 0
                    gc_counter = 0

                    if bin_chrom_start_list:  # get next startpos if there is one
                        nextpos = bin_chrom_start_list.pop(0)
                        #print(nextpos)  #TODO REMOVE WHEN DONE

                    else:  # else raise done flag and break loop
                        done_flag = 1
                        break

    gf.close()

    if outfile:
        with open(outfile, "w") as outf:
            for gc in node_gc_ratios_list:
                outf.write(str(gc) + "\n")

    return node_gc_ratios_list


# annotate graph nodes with G/C ratio as a vertex property "gc_ratio"
def annotate_nodes_with_gc(graph, gclist, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    gc_vp = g.new_vertex_property("double")

    for i in range(0, len(gclist)):
        gc_vp[g.vertex(i)] = gclist[i]

    g.vp["gc_ratio"] = gc_vp
    g.save(outfile)


# takes a graph file and annotates it with value of first principal component "ab_eigenval" for getting A/B compartments
def annotate_compartments_no_gc(graph, outfile):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    adjacency_matrix = adjacency(g, weight=g.ep.weight).toarray()  # Compute the weighted adjacency matrix
    pca = PCA(n_components=1)
    X = pca.fit_transform(adjacency_matrix)

    print(X)
    print(type(X))
    print(X.shape)
    xlist = [float(x[0]) for x in X.tolist()]
    print(xlist)
    print(len(xlist))

    ab_vp = g.new_vertex_property("double")

    for i in range(len(xlist)):
        ab_vp[g.vertex(i)] = xlist[i]

    g.vp["ab_eigenval"] = ab_vp
    g.save(outfile)


# only use if annotated with above function
def ab_cluster_report(graph, clusterfile, outfile="clusters_ab_report.tsv"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    with open(outfile, "w") as outf:
        outf.write("numV\tnumNeg\tnumPos\tposRatio\n")
        with open(clusterfile, "r") as f:
            for line in f:  # for each node in cluster
                vlist = line.split(",")
                ablist = []

                for v in vlist:  # get corresponding eigenvalues of these nodes
                    ablist.append(g.vp.ab_eigenval[v])

                numv = str(len(ablist))
                numpos = str(len([x for x in ablist if x >= 0]))
                numneg = str(len([x for x in ablist if x < 0]))
                posratio = str(float(numpos)/float(numv))

                outf.write(numv + "\t" + numneg + "\t" + numpos + "\t" + posratio + "\n")


# uses ab_eigenval and gc_ratio vertex properties to assign A/B compartment labels to nodes as vp "ab_compartment"
def determine_ab_compartments(graph, outfile, verbose=False):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    if verbose:
        # get number of nodes for use in progress bars
        vcount = 0
        for v in g.vertices():
            vcount = vcount+1

    # build list of nodes for neg and pos compartments
    neg_vlist = []
    pos_vlist = []
    if verbose:
        print("Building compartment nodelists...")
        pbar = tqdm(total=vcount)
    for v in g.vertices():
        if g.vp.ab_eigenval[v] >= 0:
            pos_vlist.append(v)
        else:
            neg_vlist.append(v)
        if verbose:
            pbar.update()

    # get gc ratios for each of those lists
    neg_gclist = []
    pos_gclist = []

    # determine which compartment is A and which is B
    if verbose:
        print("Determining compartment gc ratios...")
    for v in neg_vlist:  # get list of GC ratios for negative compartment
        neg_gclist.append(g.vp.gc_ratio[v])
    for v in pos_vlist:  # get list of GC ratios for positive compartment
        pos_gclist.append(g.vp.gc_ratio[v])

    # get the average GC content ratio for each compartment
    neg_gc_avg = sum(neg_gclist) / len(neg_gclist)
    if verbose:
        print("Negative compartment avg GC ratio is " + str(neg_gc_avg) + "...")
    pos_gc_avg = sum(pos_gclist) / len(pos_gclist)
    if verbose:
        print("Positive compartment avg GC ratio is " + str(pos_gc_avg) + "...")

    # assign "A" to whichever commpartment has the higher GC ratio, "B" to other
    if neg_gc_avg > pos_gc_avg:
        neg = "A"  # negative values are A compartment
        pos = "B"  # positive values are B compartment
        if verbose:
            print("Negative is A, Positive is B...")

    elif pos_gc_avg > neg_gc_avg:
        neg = "B"  # negative values are B compartment
        pos = "A"  # positive values are A compartment
        if verbose:
            print("Positive is A, Negative is B...")
    else:
        print("Could not determine A/B compartments. Stopping...")
        sys.exit()

    # create new vertex property "ab_compartment" and fill it with compartment labels
    ab_vp = g.new_vertex_property("string")

    if verbose:
        print("Assigning A/B compartment status to nodes...")
        pbar = tqdm(total=vcount)
    for v in g.vertices():
        if g.vp.ab_eigenval[v] >= 0:  # label positive eigenvalues with compartment determined above
            ab_vp[v] = pos
        else:  # label negative eigenvalues with compartment determined above
            ab_vp[v] = neg
        if verbose:
            pbar.update()

    # make internal
    g.vp.ab_compartment = ab_vp

    # save graph
    if verbose:
        print("Saving...")
    g.save(outfile)


# for use by annotate_ab_compartments()
# uses outputs from compartments_pca and determine_ab() to label input graph with compartments as node annotations
#TODO FINISH
def label_ab():
    x = 1


# opens tad file (.bed format) and creates a {[chrom]:[list of region bounds], ...} dictionary
def tads2dict(tadfile):

    taddict = {"chr1":[], "chr2":[], "chr3":[], "chr4":[], "chr5":[], "chr6":[], "chr7":[], "chr8":[], "chr9":[],
               "chr10":[], "chr11":[], "chr12":[], "chr13":[], "chr14":[], "chr15":[], "chr16":[], "chr17":[],
               "chr18":[], "chr19":[], "chr20":[], "chr21":[], "chr22":[], "chrX":[], "chrY":[], "chrM":[]}

    f = open(tadfile, "r")
    for line in f:
        chrom, startpos, endpos = line.split()
        binbounds = str(startpos) + "-" + str(endpos)
        taddict[chrom].append(binbounds)

    return taddict


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
    clust = OPTICS(min_samples=5, xi=0.05).fit(eigenvectors)  # do clustering
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


# performs Markov clustering on graph
def do_markov_clustering(graph, outfile, inflation=2):
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

    result = mcl.run_mcl(matrix, inflation=inflation)  # run MCL
    clusters = mcl.get_clusters(result)  # get clusters

    with open(outfile, "w") as f:
        for line in clusters:
            writeline = str(line).split("(")[-1].split(")")[0] + "\n"
            f.write(str(writeline))


# gets enrichment of GO terms for all clusters in cluster file from
def get_go_enrichments(graph, clustersfile, outfile, go_only=True):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        sys.exit()

    # read in cluster file
    with open(clustersfile, "r") as f:
        resultsDFlist = []
        for line in f:  # each line is a cluster
            splitlist = [int(x) for x in line.split(",")[:-1]]  # get list of nodes in cluster
            genelist = []
            for node in splitlist:
                genelist.extend(g.vp.genes[node])

            # for cluster in file
            gp = GProfiler(return_dataframe=True)
            if genelist:
                if go_only:
                    newDF = gp.profile(organism='hsapiens', query=genelist,
                                       sources=["GO:MF", "GO:CC", "GO:BP"], no_evidences=False)
                else:
                    newDF = gp.profile(organism='hsapiens', query=genelist, no_evidences=False)

                newDF.to_csv(outfile, mode='a')


# function for plotting GO terms / node distribution
def plot_go_terms_per_node(graphfile, binsize=5, outfile="goTermPerNodeHistogram.png"):
    g = load_graph(graphfile)  # load graph
    countlist = []  # list for storing go term count for every node
    for v in g.vertices():
        terms = g.vp.goterms[v]
        countlist.append(len(terms))  # append length of list retrieved from property map as number of terms

    # Create histogram
    plt.hist(countlist, bins=binsize, color='black', alpha=0.7, rwidth=0.85)

    # Add labels and title
    plt.xlabel('Number of GO terms per node')
    plt.ylabel('Count')

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# function for plotting genes / node distribution
def plot_genes_per_node(graphfile, binsize=5, outfile="genesPerNodeHistogram.png"):
    g = load_graph(graphfile)
    countlist = []  # list for storing go term count for every node
    for v in g.vertices():
        terms = g.vp.genes[v]
        countlist.append(len(terms))  # append length of list retrieved from property map as number of terms

    # Create histogram
    plt.hist(countlist, bins=binsize, color='black', alpha=0.7, rwidth=0.85)

    # Add labels and title
    plt.xlabel('Number of genes per node')
    plt.ylabel('Count')

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# function for plotting degree distribution
def plot_degree_distribution(graphfile, binsize=5, outfile="degreeHistogram.png"):
    g = load_graph(graphfile)
    countlist = []  # list for storing go term count for every node
    for v in g.vertices():
        countlist.append(int(v.out_degree()))  # append length of list retrieved from property map as number of terms

    # Create histogram
    plt.hist(countlist, bins=binsize, color='black', alpha=0.7, rwidth=0.85)

    # Add labels and title
    plt.xlabel('Degree per node')
    plt.ylabel('Count')

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# function for plotting distribution of edge weights
def plot_edge_weight_distribution(graphfile, binsize=100, outfile="edgeWeightHistogram.png", log=False):
    print("loading " + graphfile + " for plotting...")
    g = load_graph(graphfile)

    hist, bins, _ = plt.hist(g.ep.weight.a, bins=binsize)
    if log:
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]+1), len(bins))
    else:
        logbins = binsize

    plt.clf()
    # Create histogram
    plt.hist(g.ep.weight.a, bins=logbins, color='black')

    # Add labels and title
    plt.xlabel('Edge weight')
    plt.ylabel('Count')
    #plt.xlim(0.1, 1000)  # why does setting this result in such a strange plot?
    if log:
        plt.yscale('log')
        plt.xscale('log')  # TODO remove this line

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# function for plotting distribution of edge weights
def plot_raw_weight_distribution(graphfile, binsize=100, outfile="edgeWeightHistogram.png", log=False):
    print("loading " + graphfile + " for plotting...")
    g = load_graph(graphfile)

    hist, bins, _ = plt.hist(g.ep.weight.a, bins=binsize)
    if log:
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]+1), len(bins))
    else:
        logbins = binsize

    plt.clf()
    # Create histogram
    plt.hist(g.ep.raw_weight.a, bins=logbins, color='black')

    # Add labels and title
    plt.xlabel('Edge weight')
    plt.ylabel('Count')
    #plt.xlim(0.1, 1000)  # why does setting this result in such a strange plot?
    if log:
        plt.yscale('log')
        plt.xscale('log')  # TODO remove this line

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# function for plotting edge weights from before transformation
def plot_raw_weight_distribution(graphfile, binsize=100, outfile="rawWeightHistogram.png"):
    print("loading " + graphfile + " for plotting...")
    g = load_graph(graphfile)

    hist, bins, _ = plt.hist(g.ep.raw_weight.a, bins=binsize)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]+1), len(bins))
    plt.clf()
    # Create histogram
    plt.hist(g.ep.raw_weight.a, bins=logbins, color='black')

    # Add labels and title
    plt.xlabel('Raw edge weight')
    plt.ylabel('Count')
    plt.xlim(0.1, 1000)  # why does setting this result in such a strange plot?
    plt.yscale('log')
    plt.xscale('log')  # TODO remove this line

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# function for plotting shortest path distribution for all pairs of nodes
def plot_shortest_path_pairs_distribution(APSPmatrix, binsize=100, outfile="pairwiseShortestPathsHistogram.png", log=False):
    mat = np.loadtxt(open(APSPmatrix, "rb"), skiprows=1)[:, 1:]  # read distrn matrix w/o index col/row
    nrows, ncol = mat.shape  # get nparray shape to use as number of rows/columns

    # extract shortest paths from lower triangle of matrix
    pathweights = mat[np.tril_indices(nrows, k=-1)]  # takes values in lower triangle of mat excluding diagonal

    hist, bins, _ = plt.hist(pathweights, bins=binsize)
    if log:
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]+1), len(bins))
    else:
        logbins = binsize
    plt.clf()
    # Create histogram
    plt.hist(pathweights, bins=logbins, color='black')

    # Add labels and title
    plt.xlabel('Length(summed weight) of pairwise shortest paths')
    plt.ylabel('Count')
    if log:
        plt.xscale('log')
        plt.yscale('log')

    # Save the figure
    plt.savefig(outfile)
    plt.close()


# draws histogram of pvals in real and shuffled graphs to compare distributions (shuffled should be flat, real left skewed)
def plot_real_shuffled_pval_distributions(results_file, outfile, stepsize=0.001):
    with open(results_file, "r") as f:
        pvals = []
        shufpvals = []
        for line in f:  # extract pval and shufpval cols/rows from resultsdf (which is in a really dumb format)
            splitline = line.split(",")
            if splitline[0] == "pval":
                pvals = splitline[1:]
            elif splitline[0] == "shufpval":
                shufpvals = splitline[1:]

    pylist = [float(x) for x in pvals]  # convert to numeric
    shufpylist = [float(x) for x in shufpvals]  # convert to numeric
    xlist = list(np.arange(0, 1 + stepsize, stepsize))

    plt.clf()
    # Create histogram
    plt.hist(pylist, bins=xlist, alpha=0.5, label='real graph', color='blue', cumulative=False)
    plt.hist(shufpylist, bins=xlist, alpha=0.5, label='null graph', color='black', cumulative=False)

    # Add labels and title
    plt.xlabel('p-value', fontsize=20)
    plt.ylabel('# of GO terms in\n p-value range', fontsize=20)
    plt.xlim(xlist[0], xlist[-1])
    plt.legend(loc='upper right')
    # Save the figure
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


# draws histogram of significant pvals at decreasing thresholds to determine FDR cutoff (Panel A)
def plot_fdr_pval_histogram(results_file, outfile, stepsize=0.001, log_axis=True):
    with open(results_file, "r") as f:
        pvals = []
        shufpvals = []
        for line in f:  # extract pval and shufpval cols/rows from resultsdf (which is in a really dumb format)
            splitline = line.split(",")
            if splitline[0] == "pval":
                pvals = splitline[1:]
            elif splitline[0] == "shufpval":
                shufpvals = splitline[1:]

    pylist = [float(x) for x in pvals]  # convert to numeric
    shufpylist = [float(x) for x in shufpvals]  # convert to numeric
    xlist = list(np.arange(stepsize, 0.1 + stepsize, stepsize))

    if log_axis:
        logbins = np.logspace(np.log10(xlist[0]), np.log10(xlist[-1]), len(xlist))
    else:
        logbins = xlist

    plt.clf()
    # Create histogram
    plt.hist(pylist, bins=logbins, alpha=0.5, label='real graph', color='blue', cumulative=True)
    plt.hist(shufpylist, bins=logbins, alpha=0.5, label='null graph', color='black', cumulative=True)

    # Add labels and title
    plt.xlabel('p-value threshold', fontsize=20)
    plt.ylabel('# of GO terms with\n p-value < threshold', fontsize=20)
    plt.xlim(xlist[0], xlist[-1])  # why does setting this result in such a strange plot?
    if log_axis:
        plt.yscale('log')
        plt.xscale('log')  # TODO remove this line
    plt.legend(loc='upper left')
    # Save the figure
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


# plots correlation of TPD vs # of interchromosomal contacts for each go term or other annotation
def plot_interchromosomal_contact_tpd_correlation(resultscsv, graph, outcsv=None, outplot=None, prop="goterms"):
    # load results df to get tpds
    df = pd.read_csv(resultscsv)
    df = df.transpose()
    df.columns = df.iloc[0]
    df["goterm"] = df.index
    df = df[1:]
    df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                    "shufpval": "float", "goterm": "str"})

    # make godict to get list of nodes for each go term
    gd = make_go_dict(graph, prop=prop)  # TODO put after loading graph once I fix graph loading in all functions

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    # annotate edges with contact info
    contact_annotate(g)

    # initialize data structures
    contactsdfcol = []
    intercontactsdfcol = []
    contactratiodfcol = []


    # for each GO term, for each node, count number of interchromosomal contacts (as a function of all contacts?)
    # also fetch TPD from results df for go term while in this loop
    for index, row in df.iterrows():
        goterm = row["goterm"]
        print(goterm)
        tpd = row["tpd"]
        intercontactcounter = 0
        contactcounter = 0

        # for each node in list from godict, count intercontact edges and total edges
        for node in gd[goterm]:  # iterate through list of nodes with goterm
            v = g.vertex(node)
            for e in v.out_edges():
                if g.ep.contactType[e] == "interchromosomal":
                    intercontactcounter = intercontactcounter + 1
                contactcounter = contactcounter + 1
            contactratio = intercontactcounter / contactcounter
            print("contacts:" + str(contactcounter) + " intercontacts:" +
                  str(intercontactcounter) + " ratio:" + str(contactratio))

        # add stats to respective lists to later append to df
        contactsdfcol.append(contactcounter)
        intercontactsdfcol.append(intercontactcounter)
        contactratiodfcol.append(contactratio)

    # add to df
    df["totalContacts"] = contactsdfcol
    df["interchromosomalContacts"] = intercontactsdfcol
    df["contactRatio"] = contactratiodfcol

    # write to csv
    if outcsv is not None:
        df.to_csv(outcsv)

    # plot
    if outplot is not None:
        sns.set_theme(style="ticks")
        corplot = sns.lmplot(data=df, x="tpd", y="contactRatio",
        palette="muted", ci=None,
        height=4, scatter_kws={"s": 5, "alpha": 1})

    corplot.savefig(outplot)


# plots the number of terms passing significance threshold at each FDR value
# input should be csv from get_real_v_null_significants
def plot_passing_terms_v_fdr_fig(inputcsv, outfile="numtermsVfdrFig.png", xlim=0.25, ylim=150, logy=False):

    df = pd.read_csv(inputcsv, sep=", ")  # read input csv into data frame
    df = df.astype({"pvalThreshold": "float", "numReal": "int", "numShuf": "int", "fdr": "float"})
    newdflist = []
    # filter dataframe so there is only one row per FDR value
    currentfdr = 0
    for index, row in df.iterrows():
        if row["fdr"] != currentfdr:  # new fdr at row means we are at the new best threshold (most real hits / fdr)
            newdflist.append({"pvalThreshold": row["pvalThreshold"], "numReal": row["numReal"], "numShuf": row["numShuf"],
                              "fdr": row["fdr"]})
            currentfdr = row["fdr"]  # update to new fdr for row

    newdf = pd.DataFrame.from_records(newdflist)  # convert from list of dicts to df

    # make plot using new filtered df
    plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
    plt.scatter(newdf["fdr"], newdf["numReal"], marker='o', color='b', alpha=0.5)
    #plt.title("Number of real significant GO terms versus False Discovery Rate")
    plt.xlabel("False Discovery Rate")
    plt.xlim(0, xlim)
    if logy:
        plt.yscale("log")
    else:
        plt.ylim(0, ylim)
    plt.ylabel("Number of statistically significantly clustered GO terms")
    plt.savefig(outfile)
    plt.close()


# full run of the pipeline, from sam files to final results dataframe plus figures
# files that need to be in filedir:
# - map(HUMAN_9606_idmapping_selected.tsv)
# - gencode(gencode_pcgenes.csv)
# - go obo(go-basic.obo)
# - GothicAPSP.so
# - Ensembl2Reactome_All_Levels.txt
def gothic_full_run(runname, sam1=None, sam2=None, filedir=".", binsize=80000, vmax=10000, ncores=1,
                    step=0, endstep=12, binfile=None, saveadjlist=False, filterlevel=0):
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
        genes_go_annotate(graphname, mapfile=mymap, gencodefile=gencode, binsize=binsize,
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
        print('Global Runtime(TPDstep): ', stop - start)  # Timer Block

    # get final dataframe with "nnodes", "tpd", "pval", "shuftpd", "shufpval" for all go terms
    if step <= 11 <= endstep:
        print("Collecting final results...")
        graphname = runname + ".gt"
        distrnfile = runname + "_distributions.csv"
        resultscsvname = runname + "_resultsDF.csv"
        topgolistname = runname + "_topGOterms.csv"
        fdrthresholdsname = runname + "_FDRthresholds.csv"
        fdrplotname = runname + "_FDRpvalComparePlot.png"
        pvaldistrnplotname = runname + "_pvalHistogram.png"

        df = get_go_tpd_pvals(graphname, runname + "_realTPDs.csv", runname + "_shuffledTPDs.csv", distrnfile)
        df.to_csv(resultscsvname)

        # also get list of top go terms with FDR
        get_top_goterms(resultscsvname, outfile=topgolistname)
        get_real_v_null_significants(resultscsvname, outfile=fdrthresholdsname)

        # plots
        plot_real_shuffled_pval_distributions(resultscsvname, pvaldistrnplotname)
        plot_fdr_pval_histogram(resultscsvname, fdrplotname)

        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(EndStep): ', stop - start)  # Timer Block


def gothic_full_clustering_run(runname, sam1=None, sam2=None, filedir=".", binsize=80000, vmax=10000, ncores=1,
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
        genes_go_annotate(graphname, mapfile=mymap, gencodefile=gencode, binsize=binsize,
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

