import sys

def main():
    C01_ann_snp_file = sys.argv[1]
    E01_ann_snp_file = sys


C01_ann_snp_file = "/rd/caiya/fly/snp_sjy/C01.filter.norm.merfin.filter.afdp.ann.vcf.snp.out"
E01_ann_snp_file = "/rd/caiya/fly/snp_sjy/E01.filter.norm.merfin.filter.afdp.ann.vcf.snp.out"

C01_sp_snp_file = "/rd/caiya/fly/snp_sjy/minimap2/csp_strict1.csv.repeat.out"
E01_sp_snp_file = "/rd/caiya/fly/snp_sjy/minimap2/esp_strict1.csv.repeat.out"

out_file = path + "C01_E01.filter.norm.merfin.filter.afdp.ann.vcf.shared_sp.snp.out"
final = open(out_file, "w")

C01_sp_list = {}
with open(C01_sp_snp_file, "r") as Csp:
    for line in Csp:
        line = line.rstrip("\n")
        ll = line.split("\t")

        csp_chr = ll[0]
        csp_pos = int(ll[1])
        C01_sp_list.setdefault(csp_chr, {})[csp_pos] = 1

E01_sp_list = {}
with open(E01_sp_snp_file, "r") as Esp:
    for line in Esp:
        line = line.rstrip("\n")
        ll = line.split("\t")

        esp_chr = ll[0]
        esp_pos = int(ll[1])
        E01_sp_list.setdefault(esp_chr, {})[esp_pos] = 1

E01_total_list = {}
with open(E01_ann_snp_file, "r") as Eannnsp:
    for line in Eannnsp:
        if not line.startswith("#"):
            line = line.rstrip("\n")
            ll = line.split("\t")

            Eannsnp_chr = ll[0]
            Eannnsp_pos = int(ll[1])
            
            shared_sp_tag = "shared"
            if Eannsnp_chr in E01_sp_list:
                for esp_pos in E01_sp_list[Eannsnp_chr]:
                    if esp_pos == Eannnsp_pos:
                        shared_sp_tag = "E01_sp"
                        break
            
            E01_total_list.setdefault(Eannsnp_chr, {})[Eannnsp_pos] = line + "\t" + shared_sp_tag


C01_total_list = {}
with open(C01_ann_snp_file, "r") as Cannnsp:
    for line in Cannnsp:
        if not line.startswith("#"):
            line = line.rstrip("\n")
            ll = line.split("\t")

            Cannsnp_chr = ll[0]
            Cannnsp_pos = int(ll[1])
            
            shared_sp_tag = "shared"
            if Cannsnp_chr in C01_sp_list:
                for csp_pos in C01_sp_list[Cannsnp_chr]:
                    if csp_pos == Cannnsp_pos:
                        shared_sp_tag = "C01_sp"
                        break
            
            if shared_sp_tag != "C01_sp":
                if Cannsnp_chr in E01_total_list:
                    for Eannnsp_pos in E01_total_list[Cannsnp_chr]:
                        if Eannnsp_pos == Cannnsp_pos:
                            E01_total_list[Cannsnp_chr][Eannnsp_pos] = 0
            
            final.write(line + "\t" + shared_sp_tag + "\n")


for Eannsnp_chr in E01_total_list:
    for Eannnsp_pos in E01_total_list[Eannsnp_chr]:
        if E01_total_list[Eannsnp_chr][Eannnsp_pos] != 0:

            final.write(E01_total_list[Eannsnp_chr][Eannnsp_pos] + "\n")

final.close()

