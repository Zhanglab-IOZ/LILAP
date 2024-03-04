import sys
import pysam

def main():
    bamfile = sys.argv[1]
    outbamfile = ".".join(bamfile.split("/")[-1].split(".")[:-1]) + ".rm_inside_te_simple.bam"
    te_bed_file = "/rd/caiya/fly/snp_sjy/TE_snps/inside_te_simple_repeat/dm6.TE.s.merge.50bp.s.merge.simple_repeat.trf.s.merge.rmsk.Satellite_Simple_repeat.merge.bed"
    bed_dic = read_bed_file(te_bed_file)

    read_bam_overlap_TE(bamfile, bed_dic, outbamfile)



def snp_in_region(snp_pos, region_pos):
    spos = int(snp_pos)
    rposs = int(region_pos.split("-")[0])
    rpose = int(region_pos.split("-")[1])
    
    if spos >= rposs and spos <= rpose:
        return True
    else:
        return False

def get_cigar_length(cigar):
    tmp_num = ""
    total_length = 0
    for i in cigar:
        if i.isdigit():
            tmp_num = tmp_num + str(i)
        else:
            if i == "=" or i == "I" or i == "D" or i == "X" or i == "M":
                tmp_num_len = int(tmp_num)
                total_length += tmp_num_len
            tmp_num = ""
    return total_length

def reads_inside_te(read_pos, region_pos):    
    readposs = int(read_pos.split("-")[0])
    readpose = int(read_pos.split("-")[1])
    
    rposs = int(region_pos.split("-")[0])
    rpose = int(region_pos.split("-")[1])
    
    if readposs >= rposs and readpose <= rpose:
        return True
    else:
        return False


def read_bed_file(te_bed_file):
    bed_dic = {}
    with open(te_bed_file, "r") as bedf:
        for line in bedf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            bedpos = ll[1] + "-" + ll[2]
            if ll[0] == "chr2L" or ll[0] == "chr2R" or ll[0] == "chr3L" or ll[0] == "chr3R" or ll[0] == "chr4" or ll[0] == "chrX" or ll[0] == "chrY":
                bed_dic.setdefault(ll[0], []).append(bedpos)
    return bed_dic


def reads_inside_te(readspos, tepos):
    reads = int(readspos.split("-")[0])
    reade = int(readspos.split("-")[1])
    teposs = int(tepos.split("-")[0])
    tepose = int(tepos.split("-")[1])

    if reads >= teposs and reade <= tepose:
        return True
    else:
        return False


def read_bam_overlap_TE(bamfile, bed_dic, outbamfile):
    inbam = pysam.AlignmentFile(bamfile, "rb")
    outbam = pysam.AlignmentFile(outbamfile, "wb", template = inbam)

    for read in inbam:
        bchr = read.reference_name
        bs = read.reference_start
        be = read.reference_end
        readpos = str(bs) + "-" + str(be)

        inside_tag = "outside"
        if bchr in bed_dic:
            for bedpos in bed_dic[bchr]:
                if reads_inside_te(readpos, bedpos):
                    inside_tag = "inside"
                
        
        if inside_tag == "outside":
            outbam.write(read)
    
    inbam.close()
    outbam.close()



if __name__ == "__main__":
    main()

