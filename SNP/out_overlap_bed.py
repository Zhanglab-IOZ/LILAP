import sys

def main():
    infile = sys.argv[1]
    bedfile = sys.argv[2]
    outtag = sys.argv[3]
    outfile = ".".join(infile.split("/")[-1].split(".")[:-1]) + "." + outtag + ".out"

    bed_dic = read_bedfile(bedfile)
    read_infile(infile, bed_dic, outtag, outfile)


def read_bedfile(bedfile):
    bed_dic = {}
    with open(bedfile, "r") as bf:
        for line in bf:
            line.rstrip("\n")
            ll = line.split("\t")
            bchr = ll[0]
            bpos = ll[1] + "-" + ll[2]
            bed_dic.setdefault(bchr, []).append(bpos)
    
    return bed_dic

def snp_overlap(snp_pos, gtf_pos):
    snp_site = int(snp_pos)

    gtf_s = int(gtf_pos.split("-")[0])
    gtf_e = int(gtf_pos.split("-")[1])

    if snp_site >= gtf_s and snp_site <= gtf_e:
        return True
    else:
        return False


def read_infile(infile, bed_dic, outtag, outfile):
    final = open(outfile, "w")
    with open(infile, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            vchr = ll[0]
            snp_pos = int(ll[1])

            tag1 = "out_" + outtag
            if vchr in bed_dic:
                for bpos in bed_dic[vchr]:
                    if snp_overlap(snp_pos, bpos):
                        tag1 = "in_" + outtag
                        break
            
            final.write(line + "\t" + tag1 + "\n")
    final.close()


if __name__ == "__main__":
    main()
