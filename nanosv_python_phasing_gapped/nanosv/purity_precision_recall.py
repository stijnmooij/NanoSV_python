import re

fp = {}
tp = {}
with open("/home/cog/smooij/Desktop/true_false/tp.vcf", "r") as vcf:
    for line in vcf:
        match = re.search("PURITY_SCORE=(\d+)", line)
        for cut_off in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
            if int(match.groups(1)[0]) >= cut_off:
                if cut_off in tp:
                    tp[cut_off] += 1
                else:
                    tp[cut_off] = 1


with open("/home/cog/smooij/Desktop/true_false/fp.vcf", "r") as vcf:
    for line in vcf:
        match = re.search("PURITY_SCORE=(\d+)", line)
        for cut_off in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
            if int(match.groups(1)[0]) >= cut_off:
                if cut_off in fp:
                    fp[cut_off] += 1
                else:
                    fp[cut_off] = 1

print(tp)
print(fp)

for cut_off in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
    recall = tp[cut_off]/136
    precision = tp[cut_off]/(tp[cut_off]+fp[cut_off])
    print("\t".join([str(cut_off), str(recall), str(precision)]))