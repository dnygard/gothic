from gothic import *
import codecs

genome_file = "GRCh38_noalt_as.fa"
cases = ["1a", "1b", "1c"]
go_ab_repls_dict = {}  # dict of {"goname":{"1a":ABval}, {"1b":ABval}, {"1c":ABval}}
go_ab_array_list = []  # list of lists for later conversion to array via np.column_stack()
gonames_list = []  # list of gonames for later use as heatmap labels
go_list_file = "MEC1a1b1cSharedGOsTPD.txt"
obofile = "go-basic.obo"

# make heatmap
# read in list of go terms to plot (overlap of multiple replicates in this case)
go_list = []
with open(go_list_file, "r") as f:
    for line in f:
        go_list.append(line.strip())
print(go_list)  # TODO remove

# make goid:goname dictionary
print("creating goid:goname dictionary...")
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

# convert to array for heatmap and list of "GO:GONAME" strings to use as names
print("evoking AB array...")
for go in go_list:  # just use first dict for iterating through so order is consistent
    goname = termdict[go]  # fetch name of goterm from obo
    gonames_list.append(go + " - " + goname)  # add strings to list to be used as labels later

#np.savetxt("heatmap_vals_copy.csv", heatmap_vals_array, delimiter=",")  # save for posterity
# heatmap_file = codecs.open("heatmap_vals_copy_opp1a.csv", encoding='utf-8-sig')
heatmap_file = "heatmap_vals_copy.csv"
heatmap_vals_array = np.transpose(np.loadtxt(heatmap_file, delimiter=","))
# print(np.transpose(heatmap_vals_array))
# print(np.shape(heatmap_vals_array))

# create heatmap
vegetables = gonames_list
farmers = cases
harvest = heatmap_vals_array

fig, ax = plt.subplots()
im = ax.imshow(harvest)

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(farmers)), labels=farmers)
ax.set_yticks(np.arange(len(vegetables)), labels=vegetables)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
plt.colorbar(im)
# # Loop over data dimensions and create text annotations.
# for i in range(len(vegetables)):
#     for j in range(len(farmers)):
#         text = ax.text(j, i, harvest[i, j],
#                        ha="center", va="center", color="w")

fig.tight_layout()
plt.savefig("MEC_AB_heatmap_test.png")
plt.show()