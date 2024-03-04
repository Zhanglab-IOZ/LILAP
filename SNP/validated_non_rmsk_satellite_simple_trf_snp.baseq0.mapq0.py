import sys

def main():
    vcffile = sys.argv[1] ### "/rd/caiya/fly/snp_sjy/assembly/pav/C01_E01.asm_pav.chr.indel50bp.snp.gene.repeat.repeattype.out"
    CEtag = sys.argv[2] ###C01_E01.C01
    pileupfile = sys.argv[3]
    out_file = CEtag + "." + ".".join(vcffile.split("/")[-1].split(".")[1:-1]) + ".rm_inside_te_simple_repeat.mapq0.baseq0.out"
    non_rmsk_satellite_simple_trf_snps = read_vcf(vcffile)
    non_rmsk_satellite_simple_trf_pileup = read_pileup(pileupfile)
    validate_snp(non_rmsk_satellite_simple_trf_snps, non_rmsk_satellite_simple_trf_pileup, out_file)


def read_vcf(vcffile):
    non_rmsk_satellite_simple_trf_snps = {}
    with open(vcffile, "r") as vf:
        for line in vf:
            if not line.startswith("#"):
                line = line.rstrip("\n")
                ll = line.split("\t") ### ll[13] [14] [15]
                # if ll[14] != "Satellite" and ll[14] != "Simple_repeat" and ll[15] != "TRF_simple_repeat":
                non_rmsk_satellite_simple_trf_snps.setdefault(ll[0], {})[int(ll[1])] = line
    return non_rmsk_satellite_simple_trf_snps


def read_pileup(pileupfile):
    non_rmsk_satellite_simple_trf_pileup = {}
    with open(pileupfile, "r") as pf:
        for line in pf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            non_rmsk_satellite_simple_trf_pileup.setdefault(ll[0], {})[int(ll[1])] = ll[4]
    return non_rmsk_satellite_simple_trf_pileup

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

def validate_snp(snp_dic, pileup_dic, out_file):
    final = open(out_file, "w")
    bases = ["A", "T", "G", "C"]
    for snpchr in snp_dic:
        for nsnp in snp_dic[snpchr]:
            nsnp_line = snp_dic[snpchr][nsnp]
            ll = nsnp_line.split("\t")
            derived_base = ll[4]

            if nsnp in pileup_dic[snpchr]:
                # refallele = 0
                driallele = 0
                pileup_base = get_snp_pileup_base(pileup_dic[snpchr][nsnp])
                depth = len(pileup_base)
                for b in pileup_base:
                    b = b.upper()
                    if b in bases:
                        b = b.upper()
                        if b == derived_base:
                            driallele += 1
                
                if depth > 0 :
                    vaf = round(driallele/depth, 2)
                    if snpchr == "chrX" or snpchr == "chrY":
                        if vaf >= 0.9:
                            pileup_tag = "PASS"
                        else:
                            pileup_tag = "VAF"
                    else:
                        if vaf >= 0.3:
                            pileup_tag = "PASS"
                        else:
                            pileup_tag = "VAF"
                else:
                    vaf = 0
                    pileup_tag = "VAF"
            else:
                vaf = "*"
                pileup_tag = "error"
            
            final.write(nsnp_line + "\t" + str(vaf) + "\t" + pileup_tag + "\n")

    final.close()


if __name__ == "__main__":
    main()
