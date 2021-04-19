

with open("para_metadata2.csv", 'r') as f:
    for i in f.readlines():
        line = i.strip()
        if line:
            splitted = line.split(",")
            exon = splitted[0]
            tree = splitted[1]
            print("cp %s ." % exon)
            print("cp %s ." % tree)


