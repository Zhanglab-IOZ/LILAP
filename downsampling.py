import sys
import os
import random

def main():
    genome_size = 146638899
    infile = sys.argv[1]
    depth = int(sys.argv[2]) ####10 20 50

    fsize = os.path.getsize(infile)
    rate = (genome_size * depth) / fsize

    out_file = ".".join(infile.split("/")[-1].split(".")[:-1]) + "." + str(depth) + "x.fa"
    final = open(out_file, "w")

    downsampling(infile, rate, final)
    final.close()

def downsampling(infile, rate, final):
    with open(infile, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                c = 0
                r = random.random()
                if r <= rate:
                    c = 1
                    final.write(line)
            elif c == 1:
                final.write(line)

if __name__ == "__main__":
    main()

