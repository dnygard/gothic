from gothic import *

# for a given list of goterms, write nodelists
golist = ["GO:0004332", "GO:0005225", "GO:0015165"]
cases = ["1a"]


# make directory for go lists
directory_name = cwd = os.getcwd() + "NodeLists"
print(directory_name)
try:
    os.mkdir(directory_name)
    print(f"Directory '{directory_name}' created successfully.")
except FileExistsError:
    print(f"Directory '{directory_name}' already exists.")

for term in golist:
    for case in cases:
        graphfile = "MEC" + case + "TADgraphPt0175cut100kwindow_oneminusw.gt"
        nodelist = get_nodes_with_term(graphfile, term, cytoscape=True)
        outfile = directory_name + "/GO" + term.split(":")[1] + case + "Nodes.txt"
        with open(outfile, "w") as f:
            for node in nodelist:
                f.write(node + "\n")

