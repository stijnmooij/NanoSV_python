from utils import parse_breakpoints as breakpoint
from utils import parse_reads as read
from utils import parse_bam as bam
from utils import create_vcf as c_vcf
import NanoSV
import random
matrix = []

def print_wigard_lijst(bp):
    for segment in wigard_lijst[bp]:
        segment[1] = ";".join([str(x) for x in segment[1]])
        print("\t".join([str(i) for i in segment]))


def make_matrix(sv_id, windows):
    global matrix, wigard_lijst
    wigard_lijst = [[],[]]
    x = 0
    scores = [0,0,0,0]
    bp = -1
    for window in windows:
        matrix = []
        if x == 0:
            chr = breakpoint.structural_variants[sv_id].chr
        else:
            chr = breakpoint.structural_variants[sv_id].chr2
        bin_start = int(window[0] / NanoSV.opts_variant_bin_size)
        bin_end = int(window[1] / NanoSV.opts_variant_bin_size)
        if bin_start < 0:
            bin_start = 0
        if bin_end > c_vcf.vcf.contigs[chr][1]:
            bin_end = int(c_vcf.vcf.contigs[chr][1] / NanoSV.opts_variant_bin_size)
        for qname_clip in breakpoint.structural_variants[sv_id].ref_qname_clips[x]:
            if x == 0:
                wigard_lijst[0].append([sv_id, qname_clip, "ref", 0])
            else:
                wigard_lijst[1].append([sv_id, qname_clip, "ref", 0])
        sv_reads = get_sv_reads(sv_id, x)
        for bin in range(bin_start, bin_end+1):
            for variant_position in bam.variants[chr][bin]:
                if int(window[0]) <= variant_position <= int(window[1]):
                    matrix.append([])
                    matrix_fill_position(breakpoint.structural_variants[sv_id].ref_qname_clips[x], variant_position, sv_id, chr, bin, x)
                    matrix_fill_position(sv_reads, variant_position, sv_id, chr, bin, x)
        matrix = list(map(list, zip(*matrix)))
        if len(matrix) == 0:
            x += 1
            continue
        phasing_result= clustering(matrix, list(range(len(breakpoint.structural_variants[sv_id].ref_qname_clips[x]),len(matrix))), x)
        if sum(phasing_result[:1]) > sum(scores[:1]):
            scores = phasing_result
            scores.append(len(matrix[0]))
            bp += 1
        x += 1
        deleted = 0
        for segment in range(len(matrix)):
            if len(set(matrix[segment-deleted])) <= 1 and matrix[segment-deleted][0] == "-":
                del matrix[segment-deleted]
                deleted += 1
    print_wigard_lijst(bp)
    return scores


def get_sv_reads(sv_id, x):
    global wigard_lijst
    sv_reads = []
    for bp_id in breakpoint.structural_variants[sv_id].breakpoints:
        if x == 0:
            sv_reads.append(read.breakpoints[bp_id].segment_1['id'])
            wigard_lijst[0].append([sv_id, read.breakpoints[bp_id].segment_1['id'], "alt", 0])
        else:
            sv_reads.append(read.breakpoints[bp_id].segment_2['id'])
            wigard_lijst[1].append([sv_id, read.breakpoints[bp_id].segment_2['id'], "alt", 0])
    return sv_reads


def matrix_fill_position(sv_or_ref, variant_position, sv_id, chr, bin, x):
    global matrix
    for segment_id in sv_or_ref:
        if segment_id[2] in bam.variants[chr][bin][variant_position].segments:
            matrix[len(matrix) - 1].append(bam.variants[chr][bin][variant_position].segments[segment_id[2]][0])
        else:
            if breakpoint.structural_variants[sv_id].flag1 == 0 and x == 0:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].pos > variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')
            elif breakpoint.structural_variants[sv_id].flag1 == 16 and x == 0:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].end < variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')
            elif breakpoint.structural_variants[sv_id].flag2 == 0 and x == 1:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].end < variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')
            elif breakpoint.structural_variants[sv_id].flag2 == 16 and x == 1:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].pos > variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')


def clustering(matrix, sv_reads, bp_id):
    clustering_matrix = make_clustering_matrix(matrix)
    while len(clustering_matrix) > 2:
        keys = []
        for x in clustering_matrix:
            keys.append(x)
        highest_score = 0
        readA = 0
        readB = 0
        for i in keys:
            for key, value in clustering_matrix[i].items():
                if value >= highest_score:
                    highest_score = value
                    readA = key
                    readB = i
        if highest_score < NanoSV.opts_clustering_cutoff:
            break
        merged_name = str(readA) + "," + str(readB)
        merged_dict = {}
        for j in keys:
            if j == readA or j == readB:
                continue
            sum_of_scores = 0
            for read in [readA, readB]:
                if max(map(int, read.split(","))) >= max(map(int, j.split(","))):
                    sum_of_scores += clustering_matrix[str(read)][str(j)]
                else:
                    sum_of_scores += clustering_matrix[str(j)][str(read)]
            merged_dict[str(j)] = sum_of_scores / 2

        del_list = []
        for item, merged_value in merged_dict.items():
            if max(map(int, item.split(","))) <= max(map(int, readB.split(","))):
                continue
            clustering_matrix[str(item)][str(merged_name)] = merged_value
            del_list.append(item)
        del clustering_matrix[str(readA)]
        del clustering_matrix[str(readB)]
        for read in clustering_matrix:
            if readA in clustering_matrix[read]:
                del clustering_matrix[read][readA]
            if readB in clustering_matrix[read]:
                del clustering_matrix[read][readB]
        for item in del_list:
            del merged_dict[str(item)]
        clustering_matrix[merged_name] = merged_dict
    breakpoint_result = judge_clustering(clustering_matrix, sv_reads, len(matrix), bp_id)

    return breakpoint_result


def judge_clustering(clustering_matrix, sv_reads, total_reads, x):
    global wigard_lijst
    purity_chart = []
    phasing_chart = []
    clusters = []
    for key in clustering_matrix:
        clusters.append(key.split(","))
    longest_clusters = []
    for length_cluster in [len(clusters), len(clusters)-1]:
        long_cluster = []
        for j in range(length_cluster):
            if len(long_cluster) <= len(clusters[j]):
                long_cluster = clusters[j]
        del clusters[clusters.index(long_cluster)]
        longest_clusters.append(long_cluster)
    amounts_per_cluster = []
    clusternr = 0
    for cluster in longest_clusters:
        amounts = [0, 0]
        clusternr += 1
        for read in cluster:
            if int(read) in sv_reads:
                amounts[0] += 1
                wigard_lijst[int(x)][int(read)][2] = "alt"
            else:
                amounts[1] += 1
                wigard_lijst[int(x)][int(read)][2] = "ref"

            wigard_lijst[int(x)][int(read)][3] = clusternr
        amounts_per_cluster.append(amounts)
    pur_percentage_per_cluster = []
    phasing_percentage_per_cluster = []
    for i in range(len(amounts_per_cluster)):
        purity_percentages = []
        purity_percentages.append(amounts_per_cluster[i][0] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        purity_percentages.append(amounts_per_cluster[i][1] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        pur_percentage_per_cluster.append(purity_percentages)
        phasing_percentages = []
        phasing_percentages.append(amounts_per_cluster[i][0] / len(sv_reads) * 100)
        try:
            phasing_percentages.append(amounts_per_cluster[i][1] / (total_reads - len(sv_reads)) * 100)
        except ZeroDivisionError:
            phasing_percentages.append(0)
        phasing_percentage_per_cluster.append(phasing_percentages)


    pur_sv_score = (pur_percentage_per_cluster[0][0] - pur_percentage_per_cluster[1][0])
    pur_ref_score = (pur_percentage_per_cluster[0][1] - pur_percentage_per_cluster[1][1])
    if pur_sv_score < 0:
        pur_sv_score = pur_sv_score * -1
    if pur_ref_score < 0:
        pur_ref_score = pur_ref_score * -1
    phasing_sv_score = (phasing_percentage_per_cluster[0][0] + phasing_percentage_per_cluster[1][0])
    phasing_ref_score = (phasing_percentage_per_cluster[0][1] + phasing_percentage_per_cluster[1][1])

    purity_score = (pur_sv_score + pur_ref_score) / 2
    phasing_score = (phasing_sv_score + phasing_ref_score) / 2
    random_purities = []
    for i in range(100000):
        random_purities.append(randomise(longest_clusters, sv_reads))
    teller = 1
    for value in random_purities:
        if float(value) >= purity_score:
            teller += 1
    pvalue = teller/len(random_purities)
    purity_chart.append(int(purity_score))
    phasing_chart.append(int(phasing_score))
    return [int(purity_score), int(phasing_score), pvalue]


def randomise(longest_clusters, sv_reads):
    random_options = []
    original = [[0], longest_clusters[0],longest_clusters[1]]
    for ref in range(len(matrix) - len(sv_reads)):
        random_options.append(1)
    for alt in range(len(sv_reads)):
        random_options.append(2)
    random.shuffle(random_options)
    new_clusters = [[], [], []]
    for clusternr in range(0,3):
        if clusternr == 0:
            for zero in range(len(matrix) - (len(longest_clusters[0])+len(longest_clusters[1]))):
                new_clusters[clusternr].append(random_options[-1])
                del random_options[-1]
        else:
            for one_or_two in range(len(longest_clusters[clusternr-1])):
                new_clusters[clusternr].append(random_options[-1])
                del random_options[-1]

    pur_percentage_per_cluster = []
    amounts_per_cluster = [[new_clusters[1].count(2),new_clusters[1].count(1)], [new_clusters[2].count(2), new_clusters[2].count(1)]]
    for i in range(len(amounts_per_cluster)):
        purity_percentages = []
        purity_percentages.append(amounts_per_cluster[i][0] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        purity_percentages.append(amounts_per_cluster[i][1] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        pur_percentage_per_cluster.append(purity_percentages)
    pur_sv_score = (pur_percentage_per_cluster[0][0] - pur_percentage_per_cluster[1][0])
    pur_ref_score = (pur_percentage_per_cluster[0][1] - pur_percentage_per_cluster[1][1])
    if pur_sv_score < 0:
        pur_sv_score = pur_sv_score * -1
    if pur_ref_score < 0:
        pur_ref_score = pur_ref_score * -1

    purity_score = (pur_sv_score + pur_ref_score) / 2
    return purity_score

def make_clustering_matrix(matrix):
    clustering_matrix = {}
    for i in range(len(matrix)):
        clustering_matrix[str(i)] = {}
        for j in range(i + 1):
            mutations_in_common = 0
            if j == i:
                clustering_matrix[str(i)][str(j)] = 0
            else:
                amount_positions = len(matrix[i])
                for pos in range(len(matrix[i])):
                    if matrix[i][pos] == '-' and matrix[j][pos] == '-':
                        amount_positions -= 1
                        continue
                    if matrix[i][pos] == matrix[j][pos]:
                        mutations_in_common += 1
                if amount_positions == 0:
                    clustering_matrix[str(i)][str(j)] = 0
                else:
                    clustering_matrix[str(i)][str(j)] = mutations_in_common / amount_positions
    return clustering_matrix
