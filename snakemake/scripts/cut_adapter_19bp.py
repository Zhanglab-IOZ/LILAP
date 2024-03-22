import sys

def main():
    fasta_dic = read_infasta()
    cut_19bp_adapter(fasta_dic)

def read_infasta():
    c = 0
    fasta_dic = {}
    for line in sys.stdin:
        c += 1
        line = line.rstrip("\n")

        if c % 2 == 1:
            fasta_dic[line] = ""
        elif c % 2 == 0:
            fasta_dic[line] = fasta_dic[line] + line

    return fasta_dic

def cut_19bp_adapter(fasta_dic):
    for name in fasta_dic:
        seq = fasta_dic[name]
        subseq = seq[19:-19]
        print(name)
        print(subseq)

if __name__ == "__main__":
    main()
