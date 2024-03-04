import sys

def main():
    file_tag = sys.argv[1] ### C01_sp
    snps_file = "C01_E01.pav.min3." + file_tag + ".gt.out"
    snps_id_dic = generate_snps_id(snps_file)

    distance = ["1_9", "1_99", "1_299", "1_1000"]
    outfile = "C01_E01.pav.minreads3." + file_tag + ".gt.mnp.reads_phasing.combine.out"
    final = open(outfile, "w")
    mnp_tag_dic = open_files(distance, file_tag, snps_id_dic)
    out_max_distance(file_tag, distance, snps_id_dic, mnp_tag_dic, final)


def generate_snps_id(snps_file):
    snps_id_dic = {}
    c = 0
    with open(snps_file, "r") as sfile:
        for line in sfile:
            c += 1
            line = line.rstrip("\n")
            ll = line.split("\t")
            schr = ll[0]
            spos = int(ll[1])
            
            snp_id = "snp_" + str(c)
            snps_id_dic.setdefault(schr, {})[spos] = snp_id
    return snps_id_dic

def open_files(distance, file_tag, snps_id_dic):
    mnp_tag_dic = {}
    for d in distance:
        F = "C01_E01.pav.min3." + file_tag + ".gt.mnp_" + d + ".reads_phasing.out"
        with open(F, "r") as mnp:
            for line in mnp:
                if not line.startswith("\n"):
                    line = line.rstrip("\n")
                    ll = line.split("\t")
                    schr = ll[0]
                    spos = int(ll[1])
                    mnp_tag1 = int(ll[-1].split("-")[1])
                    
                    snp_id = snps_id_dic[schr][spos]
                    if snp_id not in mnp_tag_dic:
                        mnp_tag_dic[snp_id] = mnp_tag1
                    elif mnp_tag1 < mnp_tag_dic[snp_id]:
                        mnp_tag_dic[snp_id] = mnp_tag1
    return mnp_tag_dic


def out_max_distance(file_tag, snps_id_dic, mnp_tag_dic, final):
    max_distance_file = "C01_E01.pav.min3." + file_tag + ".gt.mnp_1_1000.reads_phasing.out"
    with open(max_distance_file, "r") as mf:
        for line in mf:
            if line.startswith("\n"):
                final.write(line)
            else:
                line = line.rstrip("\n")
                ll = line.split("\t")
                schr = ll[0]
                spos = int(ll[1])
                snp_id = snps_id_dic[schr][spos]
                mnp_tag = mnp_tag_dic[snp_id]
                outline = "\t".join(ll[:-1]) + "\t1-" + str(mnp_tag)
                final.write(outline + "\n")
    final.close()


if __name__ == "__main__":
    main()



