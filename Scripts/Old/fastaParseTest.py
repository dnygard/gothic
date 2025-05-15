# test fasta file parsing
with open("GRCh38_noalt_as.fa", "r") as f:
    line_counter = 0
    for line in f:
        if line_counter < 50:
            print(line)
            llist = [x for x in line]
            print(llist)
            
            line_counter = line_counter + 1
        else:
            break
