import sys

def main():
    infile = sys.argv[1]
    qname_pileup_file = sys.argv[2]
    d = sys.argv[3] ### 1-9 1-99 1-299 1-1000
    ds = int(d.split("-")[0])
    de = int(d.split("-")[1])

    out_file = ".".join(infile.split("/")[-1].split(".")[:-1]) + ".mnp_" + str(ds) + "_" + str(de) + ".reads_phasing.out"
    final = open(out_file, "w")

    vcf_dic = read_infile(infile)
    mnp_dic = find_mnp(vcf_dic, ds, de, d)

    snp_reads_dic = read_qname_pileup(qname_pileup_file)
    mnp_phasing_dic = same_reads_phasing_mnp(mnp_dic, snp_reads_dic)
    trim_out_mnp_dic = trim_mnp(mnp_phasing_dic, de, ds, d, vcf_dic)
    write_final(trim_out_mnp_dic, final)
    final.close()

def read_infile(infile):
    vcf_dic = {}
    with open(infile, "r") as F:
        for line in F:
            line = line.rstrip("\n")
            ll = line.split("\t")
            vchr = ll[0]
            vpos = int(ll[1])
            vcf_dic.setdefault(vchr, {})[vpos] = line
    return vcf_dic

def find_mnp(vcf_dic, ds, de, d):
    tmp = []
    out_mnp_dic = {}
    c = 0
    for vchr in sorted(vcf_dic.keys()):
        snp_pos_dic = vcf_dic[vchr]
        for pos1 in sorted(snp_pos_dic.keys()):
            pos1_list = find_mnp_1(pos1, snp_pos_dic, ds, de)

            if pos1 in tmp:
                tmp = tmp + pos1_list
                tmp = list(set(tmp))
            else:
                c += 1
                tmp_out(tmp, snp_pos_dic, out_mnp_dic, d, c)
                tmp = pos1_list
        c += 1
        tmp_out(tmp, snp_pos_dic, out_mnp_dic, d, c)
        tmp = []
    return out_mnp_dic

def find_mnp_1(pos, snp_dic, ds, de):
    pos = int(pos)
    out_pos = [pos]
    ds = int(ds)
    de1 = int(de) + 1

    for i in range(ds, de1):
        pos2 = pos + i
        if pos2 in snp_dic:
            out_pos.append(pos2)
    return out_pos

def tmp_out(tmp, snp_dic, out_dic, d, c):
    if len(tmp) > 1:
        for p in sorted(tmp):
            snp_line = snp_dic[p]
            out_dic.setdefault(c, []).append(snp_line + "\t" + d)

def get_snp_pileup_base(aaa):
    pileup_chr = ["A", "a", "T", "t", "G", "g", "C", "c", ".", ",", "*", "+", "-"]
    base = ["A", "T", "G", "C"]
    tmp = ""
    out_str = []
    snptag = 0
    aaa = aaa.upper()
    if "+" in aaa or "-" in aaa:
        snptag = 0
    else:
        snptag = 1
    if snptag == 1:
        for i in aaa:
            if i in pileup_chr:
                out_str.append(i)    
    elif snptag == 0:
        c = 0
        for i in aaa:
            if i == "+" or i == "-":
                tmp = out_str[-1] + i
                out_str = out_str[0:-1]
                numtmp = ""
                cnumtmp = 0
                c = 1
            elif c == 1:
                if i.isdigit():
                    tmp = tmp + str(i)
                    numtmp = numtmp + str(i)
                elif i in base:
                    numtmp = int(numtmp)
                    if cnumtmp < numtmp:
                        tmp = tmp + str(i)
                        cnumtmp += 1
                    else:
                        out_str.append(tmp)
                        tmp = ""
                        out_str.append(i)
                        c = 0
                        numtmp = ""
                        cnumtmp = 0
                elif i == "." or i == "," or i == "*":
                    out_str.append(tmp)
                    tmp = ""
                    out_str.append(i)
                    c = 0
            elif c == 0:
                out_str.append(i)
        if tmp != "":
            out_str.append(tmp)
    return out_str

def read_qname_pileup(qname_pileup_file):
    snp_reads_dic = {}
    with open(qname_pileup_file, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            snp_pos = ll[0] + ":" + ll[1]
            qnames_list = ll[7].split(",")
            pileup_list = get_snp_pileup_base(ll[4])
            if len(qnames_list) == len(pileup_list):
                for i in range(len(qnames_list)):
                    qname = qnames_list[i]
                    pileup_base = pileup_list[i].upper()
                    snp_reads_dic.setdefault(snp_pos, {})[qname] = pileup_base
            else:
                print("error in " + snp_pos)
    return snp_reads_dic

def if_linkage_snp_1_2(snp_pos1_cover_reads_dic, snp_pos1_support_reads_dic, snp_pos2_cover_reads_dic, snp_pos2_support_reads_dic):
    c = 0
    cover_reads_dic = {}
    for read1 in snp_pos1_cover_reads_dic:
        if read1 in snp_pos2_cover_reads_dic:
            cover_reads_dic[read1] = 1
    
    if len(cover_reads_dic) == 0:
        c += 1
    else:
        for read2 in cover_reads_dic:
            if (read2 in snp_pos1_support_reads_dic and read2 not in snp_pos2_support_reads_dic) or (read2 not in snp_pos1_support_reads_dic and read2 in snp_pos2_support_reads_dic):
                c += 1
    
    if c > 0:
        return False
    else:
        return True

def same_reads_phasing_mnp(mnp_dic, snp_reads_dic):
    mnp_dic_phase_out = {}
    c1 = 0
    tmp_dic = {}
    tmp_dic1 = {}

    for c in mnp_dic:
        mnp_list = mnp_dic[c]
        for snp_line1 in mnp_list:
            snp_ll1 = snp_line1.split("\t")
            snp_pos1 = snp_ll1[0] + ":" + snp_ll1[1]
            if snp_pos1 not in tmp_dic1:
                der_base1 = snp_ll1[4]
                snp_pos1_cover_reads_dic = snp_reads_dic[snp_pos1]
                snp_pos1_support_reads_dic = {}
                for snp_read in snp_pos1_cover_reads_dic:
                    if snp_pos1_cover_reads_dic[snp_read] == der_base1:
                        snp_pos1_support_reads_dic[snp_read] = 1
                
                for snp_line2 in mnp_list:
                    snp_ll2 = snp_line2.split("\t")
                    snp_pos2 = snp_ll2[0] + ":" + snp_ll2[1]
                    if snp_pos2 != snp_pos1 and snp_pos2 not in tmp_dic1:
                        der_base2 = snp_ll2[4]
                        snp_pos2_cover_reads_dic = snp_reads_dic[snp_pos2]
                        snp_pos2_support_reads_dic = {}
                        for snp_read in snp_pos2_cover_reads_dic:
                            if snp_pos2_cover_reads_dic[snp_read] == der_base2:
                                snp_pos2_support_reads_dic[snp_read] = 1
                        
                        if if_linkage_snp_1_2(snp_pos1_cover_reads_dic, snp_pos1_support_reads_dic, snp_pos2_cover_reads_dic, snp_pos2_support_reads_dic):
                            tmp_dic1[snp_pos1] = snp_line1
                            tmp_dic1[snp_pos2] = snp_line2
                            tmp_dic.setdefault(snp_pos1, []).append(snp_pos2)
    
    return tmp_dic

def trim_mnp(mnp_phasing_dic, de, ds, d, vcf_dic):
    out_mnp_dic = {}
    c = 0
    for mnp_snp_1 in mnp_phasing_dic:
        tmp_snp_pos_dic = {}
        mnp_chr = mnp_snp_1.split(":")[0]
        snp_pos_dic = vcf_dic[mnp_chr]

        mnp_snp_pos1 = int(mnp_snp_1.split(":")[1])
        tmp_snp_pos_dic[mnp_snp_pos1] = 1
        for mnp_snp_2 in mnp_phasing_dic[mnp_snp_1]:
            mnp_snp_pos2 = int(mnp_snp_2.split(":")[1])
            tmp_snp_pos_dic[mnp_snp_pos2] = 1

        tmp = []
        for pos1 in sorted(tmp_snp_pos_dic.keys()):
            pos1_list = find_mnp_1(pos1, tmp_snp_pos_dic, ds, de)
            if pos1 in tmp:
                tmp = tmp + pos1_list
                tmp = list(set(tmp))
            else:
                c += 1
                tmp_out(tmp, snp_pos_dic, out_mnp_dic, d, c)
                tmp = pos1_list
        c += 1
        tmp_out(tmp, snp_pos_dic, out_mnp_dic, d, c)
        tmp = []
    
    return out_mnp_dic

def write_final(trim_out_mnp_dic, final):
    for c in trim_out_mnp_dic:
        for snp_line in trim_out_mnp_dic[c]:
            final.write(snp_line + "\n")
        final.write("\n")

if __name__ == "__main__":
    main()

