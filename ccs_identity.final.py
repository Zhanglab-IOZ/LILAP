import sys
import math

def main():
    samfile = sys.argv[1]
    tag1 = sys.argv[2] ###C01 or E01
    tag2 = sys.argv[3] ###asm maskedasm dm6 maskeddm6
    zmw_file = tag1 + ".zmw.count.out"

    out1 = zmw_file + "." + tag2 + ".mapq_more_0.out1"
    final1 = open(out1, "w")


    zmw_dic = read_zmw_count_file(zmw_file)
    read_samfile_get_qv(samfile, zmw_dic, final1)
    final1.close()


def read_zmw_count_file(zmw_file):
    zmw_dic = {}
    with open(zmw_file, "r") as zz:
        for line in zz:
            line = line.rstrip("\n")
            ll = line.split("\t")
            zmw_dic[ll[0]] = line
    return zmw_dic


def read_samfile_get_qv(samfile, zmw_dic, final1):
    with open(samfile, "r") as ss:
        for line in ss:
            line = line.rstrip("\n")
            ll = line.split("\t")
            name = ll[0].rstrip("/ccs")

            if name in zmw_dic and int(ll[4]) > 0:
                cigar = ll[5]
                if cigar == "*":
                    continue
                else:
                    rr = get_cigar(cigar)
                    final1.write(zmw_dic[name] + "\t" + str(rr.split("-")[0]) + "\t" + str(rr.split("-")[1]) + "\n")


def get_cigar(cigar):
    tmp_num = ""
    total_length = 0
    identical_length = 0
    for i in cigar:
        if i.isdigit():
            tmp_num = tmp_num + str(i)
        else:
            if i == "=" or i == "I" or i == "D" or i == "X":
                tmp_num_len = int(tmp_num)    
                total_length += tmp_num_len
            if i == "=":
                identical_length += tmp_num_len
            tmp_num = ""

    rr = str(identical_length) + "-" + str(total_length)
    return rr

def get_mean(ll):
    mm = sum(ll)/len(ll)
    return mm

def get_median(ll):
    ll.sort()
    length = len(ll)
    if length == 0:
        return 0
    if length % 2 == 1:
        md = ll[length // 2]
    else:
        md = (ll[length // 2] + ll[(length // 2) -1]) / 2
    return md

def get_qv(a, b):
    rr_1 = 1/(1+b)
    rr = float(int(a)/int(b))
    ee = float(1-rr)
    if a == b:
        qv_1 = (-10) * math.log10(rr_1)
        return qv_1
    else:
        qv_0 = (-10) * math.log10(ee)
        qv_1 = (-10) * math.log10(rr_1) 
        if qv_0 <= qv_1:
            return qv_0
        else:
            return qv_1

if __name__ == "__main__":
    main()

