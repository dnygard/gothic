from gothic import *

for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
    print(case)
    run = case + "FixedBins500kb"
    infile = run + "_COO.tsv"
    outfile = run + "_NORMEDadjlist.tsv"

    ice_balance(infile, outfile)
