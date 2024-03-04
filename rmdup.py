import sys

def main():
    samfile = sys.argv[1]
    
    out1 = ".".join(samfile.split("/")[-1].split(".")[:-1]) + "rmdup.out"
    final1 = open(out1, "w")

    read_samfiles(samfile, final1)
    final1.close()

def read_samfiles(samfile, final1):
    sam_out = {}
    with open(samfile, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                final1.write(line)
            else:
                line = line.rstrip("\n")
                ll = line.split("\t")
                if ll[0] in sam_out:
                    sam_out[ll[0]] = compare_line(line, sam_out[ll[0]])
                else:
                    sam_out[ll[0]] = line
    for i in sam_out:
        final1.write(sam_out[i] + "\n")

def compare_line(la, lb):
    lla = la.split("\t")
    llb = lb.split("\t")
    if int(lla[4]) == int(llb[4]):
        if len(lla[5]) <= len(llb[5]):
            return la
        else:
            return lb
    elif int(lla[4]) > int(llb[4]):
        return la
    elif int(lla[4]) < int(llb[4]):
        return lb

if __name__ == "__main__":
    main()



