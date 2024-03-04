file = "D:\AA\SV/non-B/complex.fas.all.out"

with open(file, 'r') as f:
    content = f.readlines()
d = {}
for line in content:
    line = line.split('\t')
    if line[0] in list(d.keys()):
        d[line[0]].append(line[2])
    else:
        d[line[0]] = [line[2]]
with open(file + '.sum', 'a') as f:
    f.write("id\tA_Phased_Repeat\tG_Quadruplex_Motif\tDirect_Repeat\tInverted_Repeat\tMirror_Repeat\tShort_Tandem_Repeat\tZ_DNA_Motif\tAll\n")
for key, value in d.items():
    s = key
    l = (key.split(':')[1]).split('-')
    l = int(l[1])-int(l[0]) +2000
    n1 = (value.count('A_Phased_Repeat')) * 2000 / l
    n2 = (value.count('G_Quadruplex_Motif')) * 2000 / l
    n3 = (value.count('Direct_Repeat')) * 2000 / l
    n4 = (value.count('Inverted_Repeat')) * 2000 / l
    n5 = (value.count('Mirror_Repeat')) * 2000 / l
    n6 = (value.count('Short_Tandem_Repeat')) * 2000 / l
    n7 = (value.count('Z_DNA_Motif')) * 2000 / l
    n8 = n1+n2+n3+n4+n5+n6+n7
    with open(file + '.sum', 'a') as f:
        f.write(key+'\t'+str(n1)+'\t'+str(n2)+'\t'+str(n3)+'\t'+str(n4)+'\t'+str(n5)+'\t'+str(n6)+'\t'+str(n7)+'\t'+str(n8)+'\n')