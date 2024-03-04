import sys

def main():
    infile = sys.argv[1]
    out_file = ".".join(infile.split("/")[-1].split(".")[:-1]) + ".vaf.out"
    final = open(out_file, "w")
    pileupfile = sys.argv[3]

    vcf_dic = read_vcf(infile)
    pileup_dic = read_pileup(pileupfile)
    calculate_vaf1(vcf_dic, pileup_dic, final)


def read_vcf(vcffile):
    vcf_dic = {}
    with open(vcffile, "r") as vf:
        for line in vf:
            if not line.startswith("#"):
                line = line.rstrip("\n")
                ll = line.split("\t") ### ll[13] [14] [15]
                # if ll[14] != "Satellite" and ll[14] != "Simple_repeat" and ll[15] != "TRF_simple_repeat":
                vcf_dic.setdefault(ll[0], {})[int(ll[1])] = line
    return vcf_dic

def read_pileup(pileupfile):
    pileup_dic = {}
    with open(pileupfile, "r") as pf:
        for line in pf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            pileup_dic.setdefault(ll[0], {})[int(ll[1])] = ll[4]
    return pileup_dic

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

def calculate_vaf0(deri, pieup_str):
    pileup_list = get_snp_pileup_base(pieup_str)
    total_count = 0
    deri_count = 0
    for i in pileup_list:
        i = i.upper()
        total_count += 1
        if i == deri:
            deri_count += 1
    if total_count != 0:
        vaf = round(deri_count/total_count, 2)
    else:
        vaf = 0
    out_str = str(vaf) + "\t" + str(total_count) + "\t" + str(deri_count)
    return out_str

def calculate_vaf1(vcf_dic, pileup_dic, final):
    for schr in vcf_dic:
        for spos in vcf_dic[schr]:
            if spos in pileup_dic[schr]:
                pileup_str = pileup_dic[schr][spos]
                snpline = vcf_dic[schr][spos]
                ll = snpline.split("\t")
                dri = ll[4]
                vaf_out_str = calculate_vaf0(dri, pileup_str)
            else:
                vaf_out_str = "*\t*\t*"
            final.write(snpline + "\t" + vaf_out_str + "\n")
    final.close()


if __name__ == "__main__":
    main()


