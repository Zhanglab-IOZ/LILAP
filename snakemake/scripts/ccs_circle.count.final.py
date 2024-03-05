import sys

def main():
    subreas_fasta_file = sys.argv[1] ####absolute path to subreads fasta file
    ccs_fasta_file = sys.argv[2] ###absolute path to ccs reads fasta file
    prefix = sys.argv[3] ##C01 E01

    out1 = prefix + ".zmw.count.out"
    out2 = prefix + ".ccs.circle.count.out"
    final1 = open(out1, "w")
    final2 = open(out2, "w")

    ccs_reads_dic = read_ccs_file(ccs_fasta_file)
    zmw_count_dic = read_subreads_file(subreas_fasta_file, ccs_reads_dic)
    total_count = write_zmw_count1(zmw_count_dic, final1)
    write_zmw_count2(total_count, final2)
    final1.close()
    final2.close()


def read_ccs_file(ccs_fasta_file):
    ccs_reads_dic = {}
    with open(ccs_fasta_file, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                line = line.rstrip("\n")
                ccs_reads_dic[line.split("/")[1]] = 1
    return ccs_reads_dic

def read_subreads_file(subreas_fasta_file, ccs_reads_dic):
    zmw_count_dic = {}
    with open(subreas_fasta_file, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                line = line.rstrip("\n")
                line = line.lstrip(">")
                if line.split("/")[1] in ccs_reads_dic:
                    zmw = "/".join(line.split("/")[0:2])
                    if zmw in zmw_count_dic:
                        zmw_count_dic[zmw] += 1
                    else:
                        zmw_count_dic[zmw] = 1
    return zmw_count_dic

def write_zmw_count1(zmw_count_dic, final1):
    total_count = {}
    for zmw in zmw_count_dic:
        c = zmw_count_dic[zmw]
    final1.write(zmw+"\t"+str(c)+"\n")
    total_count[c] += 1
    return total_count

def write_zmw_count2(total_count, final2):
    max_circle = max(total_count.keys()) + 1
    for i in range(max_circle):
        if i in total_count:
            final2.write(str(i) + "\t" + str(total_count[i]) + "\n")
        elif i != 0:
            final2.write(str(i) + "\t0\n")


if __name__ == "__main__":
    main()



