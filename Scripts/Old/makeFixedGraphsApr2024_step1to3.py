from gothic import *

for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
    print(case)
    samone = case + "_1.sam"
    samtwo = case + "_2.sam"
    run = case + "FixedBins500kb"
    filefolder = "/Users/dallasnygard/UOdrive/GothicFiles"

    gothic_full_run(run, samone, samtwo, binsize=500000, vmax=10000, 
saveadjlist=True, filedir=filefolder, endstep=3)
    
