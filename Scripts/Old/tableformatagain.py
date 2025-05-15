# reformat TPD table because it has erroneous line breaks for some reason

file = "MEC1cTADgraphPt0175cut100kwindow_TPD_SupplementalTable.csv"
writefile = "MEC1cTADgraphPt0175cut100kwindow_TPD_SupplementalTable.csv"

writelist = []
lineflag = 0
with open(file, "r") as f:
    for line in f:
        if line.startswith("GOterm"):
            writelist.append(line)
            continue
        if line.startswith("GO:") and len(line.split(",")) > 6:
            writelist.append(line)
            continue
        if lineflag:
            writelist.append(firststring + line)
            lineflag = 0
        else:
            firststring = line.strip()
            lineflag = 1

with open(writefile, "w") as f:
    for line in writelist:
        f.write(line)
