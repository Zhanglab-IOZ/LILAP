import sys
import pysam

def main():
    bamfile = sys.argv[1]
    outfile = ".".join(bamfile.split("/")[-1].split(".")[:-1]) + ".Fsa.bam"
    
    reads_counts = read_name_count(bamfile)
    read_name_count2(bamfile, reads_counts, outfile)


def read_name_count(bamfile):
    infile = pysam.AlignmentFile(bamfile, "rb")
    count = {}
    for read in infile:
        if read.qname in count:
            count[read.qname] += 1
        else:
            count[read.qname] = 1
    infile.close()
    return count

def read_name_count2(bamfile, reads_counts, outfile):
    infile = pysam.AlignmentFile(bamfile, "rb")
    final = pysam.AlignmentFile(outfile, "wb", template = infile)
    for read in infile:
        if reads_counts[read.qname] == 1:
            final.write(read)
    infile.close()
    final.close()





if __name__ == "__main__":
    main()
    
