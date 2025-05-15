# take a list of GO ids, add GO terms to table

obofile = "go-basic.obo"
resultsfile = "MECtadgraphPt0175cut100kwindow_topGOterms.csv"
outfile = "MECtadgraphPt0175cut100kwindow_topGOtermsWnames.csv"
graphfile = "MECtadgraphPt0175cut100kwindow.gt"

# make goid:goname dictionary
with open(obofile, "r") as gof:
    termdict = {}  # dictionary for mapping go ids to terms (key=ID, value=term)
    ids = []  # list to append ids and alt ids to
    for line in gof:
        splitline = line.split(": ")
        if splitline[0] == "name":
            goname = splitline[1].strip()
        if splitline[0] == "id":
            ids.append(splitline[1].strip())
        if splitline[0] == "altid":
            ids.append(splitline[1].strip())
        if splitline[0] == "def":
            for goid in ids:
                termdict[goid] = goname
            ids = []  # empty list of ids

# for key in termdict.keys():
#     print(key + "->" + termdict[key])
print(len(termdict.keys()))
print(termdict.keys())

with open(resultsfile, "r") as rf:
    with open(outfile, "w") as outf:
        outf.write("goterm, pval, FDRatThisCutoff, goname, PercentCompartmentA\n")
        for line in rf:
            if line.startswith(","):
                continue
            splitline = line.split(",")
            goterm = splitline[1].strip()
            pval = splitline[2]
            fdr = splitline[3].strip()
            if goterm in termdict.keys():
                goname = termdict[goterm]
            else:
                goname = "Unknown"

            # get A/B compartment ratio of all nodes in cluster
            vlist = get_nodes_with_term(graphfile, goterm)
            g = load_graph(graphfile)
            ablist = [g.vp.ab_compartment[i] for i in vlist]
            abratio = sum(i == "A" for i in ablist)/len(ablist)

            outf.write(goterm + ", " + pval + ", " + fdr + ", " + goname + ", " + abratio + "\n")


