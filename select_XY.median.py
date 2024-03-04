import sys
import numpy as np
from decimal import Decimal

def main():
    tag = sys.argv[1] ###C01 or E01 or dm6
    level = sys.argv[2] ###contigs or scaffolds dm6
    contig_fasta_file = sys.argv[3] ### assembly file
    male_depth = "C01_6.ccs.male." + tag + "." + level + ".depth"
    female_depth = "C01_6.ccs.female." + tag + "." + level + ".depth"

    out1 = tag + "." + level + "." + "X_candidates.median.out"
    out2 = tag + "." + level + "." + "Y_candidates.median.out"
    out3 = tag + "." + level + "." + "abnormal_contigs.median.out"
    final1 = open(out1, "w")
    final2 = open(out2, "w")
    final3 = open(out3, "w")

    out4 = tag + "." + level + "." + "X_candidates.median.fasta"
    out5 = tag + "." + level + "." + "Y_candidates.median.fasta"
    final4 = open(out4, "w")
    final5 = open(out5, "w")


    male_relative_depth_dic = read_get_relative_contig_depth(male_depth)
    female_relative_depth_dic = read_get_relative_contig_depth(female_depth)
    X_Y_candidates_dic = find_X_Y_related_contigs(male_relative_depth_dic, female_relative_depth_dic, final1, final2, final3)
    get_X_Y_candidates_fasta(contig_fasta_file, X_Y_candidates_dic, final4, final5)

    final1.close()
    final2.close()
    final3.close()
    final4.close()
    final5.close()


def read_get_relative_contig_depth(male_female_depth):
    chr_depth_dic = {}
    median_list = []
    with open(male_female_depth, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            chr_depth_dic.setdefault(ll[0], []).append(int(ll[2]))

            dd = int(ll[2])
            median_list.append(dd)

        med_depth = np.median(median_list)
        relative_depth_dic = {}
        for contig in sorted(chr_depth_dic):
            median_contig = np.median(chr_depth_dic[contig])
            rdepth = Decimal((median_contig/med_depth) * 2).quantize(Decimal("1"), rounding = "ROUND_HALF_UP")
            relative_depth_dic[contig] = rdepth
        
    return relative_depth_dic    

def find_X_Y_related_contigs(male_relative_depth_dic, female_relative_depth_dic, final1, final2, final3):
    X_Y_candidates_dic = {}
    for contig in male_relative_depth_dic:
        if contig not in female_relative_depth_dic:
            female_relative_depth_dic[contig] = 0

        out_str = contig + "\t" + str(male_relative_depth_dic[contig]) + "\t" + str(female_relative_depth_dic[contig])
        if int(male_relative_depth_dic[contig]) == 1:
            if int(female_relative_depth_dic[contig]) == 0:
                final2.write(out_str + "\n")
                X_Y_candidates_dic[contig] = "Y"
            elif int(female_relative_depth_dic[contig]) == 2:
                final1.write(out_str + "\n")
                X_Y_candidates_dic[contig] = "X"
            else:
                final3.write(out_str + "\tmaybe_Y_condidates\n")
        elif int(male_relative_depth_dic[contig]) > 2 or int(male_relative_depth_dic[contig]) == 0:
            final3.write(out_str + "\n")
    
    return X_Y_candidates_dic



def get_X_Y_candidates_fasta(contig_fasta_file, X_Y_candidates_dic, final4, final5):
    c = 0
    with open(contig_fasta_file, "r") as fasta:
        for line in fasta:
            line = line.rstrip("\n")
            if line.startswith(">"):
                c = 0
                seqname = line.lstrip(">")
                if seqname in X_Y_candidates_dic:
                    if X_Y_candidates_dic[seqname] == "X":
                        c = 1
                    elif X_Y_candidates_dic[seqname] == "Y":
                        c = 2
                        
            if c == 1:
                final4.write(line+"\n")
            elif c == 2:
                final5.write(line+"\n")


if __name__ == "__main__":
    main()

