import re
from utils import parse_breakpoints as breakpoint
import NanoSV
from random import shuffle


def make_total_matrix(svs_and_segments, readnames_sv_db):
    verified_sv_ids = [10224, 25245, 10282, 10304, 10228, 10305, 10237, 10311, 29116, 25790, 25259, 29118, 4781, 10509, 10261, 10511, 25283, 10266, 25801, 29119, 10345, 8020, 10367, 25810, 10218, 31987, 32428, 32523, 32511, 32000, 32020, 11087, 25302, 12604, 31734, 32076, 32086, 25355, 31854, 24474]
    verified_sv_ids_old = [24801, 9939, 9996, 10564, 9943, 10566, 9952, 10573, 28626, 25343, 24815, 28628, 4634, 10332,
                       9975, 10338, 24839, 9980, 25357, 28629, 10608, 7771, 10630, 25366, 9933, 32030, 32521, 32721,
                       32709, 32043, 32065, 11135, 24860, 12601, 31239, 32125, 32136, 24911, 31709, 24055]
    verified_sv_id_random = [2435, 16122, 27001, 2196, 25590, 26749, 26992, 26024, 27884, 2902, 17207, 26927, 4873,
                             22200, 28529, 24999, 30691, 25695, 14488, 23223, 6181, 19045, 9570, 25705, 18732, 23653,
                             26307, 25610, 27208, 32238, 32392, 28134, 8975, 10608, 12881, 21694, 23258, 2341, 9579]
    verified_sv_id = [10224]
    #sv 10224= 2:167032588
    print("make_cluster_matrix aangeroepen")
    for sv_id in svs_and_segments:
        sv_result = []
        bar_chart_score = []
        if sv_id not in verified_sv_id or breakpoint.structural_variants[sv_id].format['GT'] == '0/0':
            continue
        # if sv_id == 10224:
        #     print("Segmenten", svs_and_segments[sv_id])
        print("\n\n")
        breakpoint_list_id = 0
        for breakpoint_segments in svs_and_segments[sv_id]:
            positions_set = set()
            # print("booleans", breakpoint.structural_variants[sv_id].flag1 == 0, breakpoint.structural_variants[sv_id].flag2 == 16)
            # TT orientation
            if breakpoint.structural_variants[sv_id].flag1 == 16 and breakpoint.structural_variants[sv_id].flag2 == 0:
                for segment in breakpoint_segments:
                    for position in segment.variations:
                        if breakpoint_list_id == 0 and position > breakpoint.structural_variants[
                            sv_id].pos - NanoSV.opts_phasing_window and position < breakpoint.structural_variants[sv_id].pos:
                            positions_set.add(position)
                        elif breakpoint_list_id == 1 and position > breakpoint.structural_variants[sv_id].info[
                            'END'] - NanoSV.opts_phasing_window and position <= breakpoint.structural_variants[sv_id].info['END']:
                            positions_set.add(position)
            # HH orientation
            elif breakpoint.structural_variants[sv_id].flag1 == 0 and breakpoint.structural_variants[sv_id].flag2 == 16:
                for segment in breakpoint_segments:
                    for position in segment.variations:
                        if breakpoint_list_id == 0 and breakpoint.structural_variants[sv_id].pos <= position < breakpoint.structural_variants[sv_id].pos + NanoSV.opts_phasing_window:
                            positions_set.add(position)
                        elif breakpoint_list_id == 1 and breakpoint.structural_variants[sv_id].info['END'] <= position < breakpoint.structural_variants[sv_id].info['END'] + NanoSV.opts_phasing_window:
                            positions_set.add(position)
            # TH orientation
            elif breakpoint.structural_variants[sv_id].flag1 == 16 and breakpoint.structural_variants[sv_id].flag2 == 16:
                for segment in breakpoint_segments:
                    for position in segment.variations:
                        if breakpoint_list_id == 0 and position > breakpoint.structural_variants[
                            sv_id].pos - NanoSV.opts_phasing_window and position <= breakpoint.structural_variants[sv_id].pos:
                            positions_set.add(position)
                        elif breakpoint_list_id == 1 and position >= breakpoint.structural_variants[sv_id].info[
                            'END'] and position < breakpoint.structural_variants[sv_id].info['END'] + NanoSV.opts_phasing_window:
                            positions_set.add(position)
            # HT orientation
            elif breakpoint.structural_variants[sv_id].flag1 == 0 and breakpoint.structural_variants[sv_id].flag2 == 0:
                for segment in breakpoint_segments:
                    for position in segment.variations:
                        if breakpoint_list_id == 0 and position >= breakpoint.structural_variants[sv_id].pos and position < \
                                breakpoint.structural_variants[sv_id].pos + NanoSV.opts_phasing_window:
                            positions_set.add(position)
                        elif breakpoint_list_id == 1 and position > breakpoint.structural_variants[sv_id].info[
                            'END'] - NanoSV.opts_phasing_window and position <= breakpoint.structural_variants[sv_id].info['END']:
                            positions_set.add(position)

            positions = sorted(positions_set)
            # print("SV_id", sv_id)
            # if sv_id == 10224:
                # print("positions", positions)
            matrix = []
            readnames = []
            sv_reads = []
            readnames_db = {}
            for segment in svs_and_segments[sv_id][breakpoint_list_id]:
                segment_matrix = []
                readnames.append(segment.qname)
                # geen var: 167032769
                # print("positions", positions)
                if len(matrix) > 14:
                    print(segment.variations, "-1 check")
                for position in positions:
                    try:
                        segment_matrix.append(segment.variations[position][0])
                    except KeyError:
                        # print(position, segment.pos, segment.end, "pos, start, end", segment.qname)
                        if position >= segment.pos and position <= segment.end:
                            segment_matrix.append(5)
                        else:
                            segment_matrix.append(None)
                if segment.qname in readnames_sv_db[sv_id]:
                    readnames_db[len(matrix)] = segment.qname
                    sv_reads.append(len(matrix))
                matrix.append(segment_matrix)
                # if sv_id == 10224:
                #     print("10224 matrix", segment_matrix)
            to_be_deleted = []
            for pos in range(len(positions)):
                values_of_position_highq = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
                values_of_position_lowq = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
                total_reads_on_position = 0
                for i in range(len(matrix)):
                    try:
                        if matrix[i][pos] == 5:
                            values_of_position_highq[5] += 1
                            total_reads_on_position += 1
                        elif svs_and_segments[sv_id][breakpoint_list_id][i].variations[positions[pos]] == [0]:
                            values_of_position_highq[0] += 1
                            total_reads_on_position += 1
                        elif matrix[i][pos] is not None and \
                                svs_and_segments[sv_id][breakpoint_list_id][i].variations[positions[pos]][1] >= \
                                NanoSV.opts_min_base_qual_ph:
                            values_of_position_highq[matrix[i][pos]] += 1
                            total_reads_on_position += 1
                        elif svs_and_segments[sv_id][breakpoint_list_id][i].variations[positions[pos]][1] < \
                                NanoSV.opts_min_base_qual_ph:
                            values_of_position_lowq[matrix[i][pos]] += 1
                            total_reads_on_position += 1
                    except (KeyError, IndexError) as e:
                        continue
                true_calls = 0
                if values_of_position_highq[0] > (NanoSV.opts_max_deletions * total_reads_on_position):
                    true_calls -= 5
                for key, value in values_of_position_highq.items():
                    if key == 0:
                        continue
                    if key != 5 and value > NanoSV.opts_min_occurences_of_var * total_reads_on_position and value + values_of_position_lowq[key] < int(NanoSV.opts_max_occurences_of_var * total_reads_on_position):
                        true_calls += 1
                        # print(values_of_position_highq)
                        # print(values_of_position_lowq)
                largest = 0
                second_largest = 0
                print(values_of_position_highq, "highq")
                for key, value in values_of_position_highq.items():
                    if value > largest:
                        second_largest = largest
                        largest = key
                for i in range(len(matrix)):
                    print("matrix [i]", matrix[i])
                    if matrix[i][pos] != largest and matrix[i][pos] != second_largest:
                        matrix[i][pos] = -1
                if true_calls < 1:
                    to_be_deleted.append(positions[pos])
            print("sv_id", sv_id, "aantal reads in sv", len(svs_and_segments[sv_id][breakpoint_list_id]))
            if breakpoint_list_id == 0:
                print("chr", breakpoint.structural_variants[sv_id].chr, "pos", breakpoint.structural_variants[sv_id].pos)
            else:
                print("chr", breakpoint.structural_variants[sv_id].chr2, "pos", breakpoint.structural_variants[sv_id].info['END'])
            if len(positions) > len(to_be_deleted) and len(sv_reads) >= 0.3 * (len(matrix) - len(sv_reads)):
                breakpoint_result = filter_matrix(matrix, positions, to_be_deleted, sv_reads, bar_chart_score)
                sv_result.append([breakpoint_result[0], breakpoint_result[1]])
            elif len(positions) == len(to_be_deleted):
                print("geen informatieve SNP-posities")
                sv_result.append([0, 0])
                continue
            elif len(sv_reads) < 0.3 * (len(matrix) - len(sv_reads)):
                print("te weinig sv reads")
                print("SV reads:", sv_reads)
                sv_result.append([0, 0])
                continue
            breakpoint_list_id += 1
            print("SV reads:", sv_reads)
            print(sv_result)
        try:
            if sum(sv_result[0]) >= sum(sv_result[1]):
                bar_chart_score.append(sv_result[0])
            else:
                bar_chart_score.append(sv_result[1])
        except IndexError:
            bar_chart_score.append(sv_result[0])
        print("")
    header_string = ""
    purity_string = ""
    phasing_string = ""
    for sv in range(len(bar_chart_score)):
        header_string += str(sv) + "\t"
        purity_string += str(bar_chart_score[sv][0]) + "\t"
        phasing_string += str(bar_chart_score[sv][1]) + "\t"

    file = open("r_file.dat", "w")
    file.write(header_string + "\n")
    file.write(purity_string + "\n")
    file.write(phasing_string + "\n")
    file.close()


def print_base(hooglaag, nummer):
    if nummer == 0:
        print(hooglaag, "del")
    elif nummer == 1:
        print(hooglaag, "A")
    elif nummer == 2:
        print(hooglaag, "C")
    elif nummer == 3:
        print(hooglaag, "G")
    elif nummer == 4:
        print(hooglaag, "T")
    elif nummer == 5:
        print(hooglaag, "ref")


def filter_matrix(matrix, positions, to_be_deleted, sv_reads, bar_chart_score):
    for pos in to_be_deleted:
        index = positions.index(pos)
        for segment in matrix:
            del segment[index]
        del positions[index]
    print(positions, "overgebleven posities")
    print(matrix)
    breakpoint_result = clustering(matrix, sv_reads, bar_chart_score)
    return breakpoint_result


def clustering(matrix, sv_reads, bar_chart_score):
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

    breakpoint_result = judge_clustering(clustering_matrix, sv_reads, len(matrix), bar_chart_score)

    print("Clustering done")
    return breakpoint_result


def judge_clustering(clustering_matrix, sv_reads, total_reads, bar_chart_score):
    purity_chart = []
    phasing_chart = []
    sv_reads = [[i] for i in range(total_reads)]
    shuffle(sv_reads)
    sv_reads = sv_reads[:int(len(sv_reads)/2)]
    bad = False
    clusters = []
    for key in clustering_matrix:
        print(key)
        clusters.append(key.split(","))
        if ',' not in key:
            print("slechte clustering")
            bad = True
    longest_clusters = []
    for length_cluster in [len(clusters), len(clusters)-1]:
        long_cluster = []
        for j in range(length_cluster):
            if len(long_cluster) <= len(clusters[j]):
                long_cluster = clusters[j]
        del clusters[clusters.index(long_cluster)]
        longest_clusters.append(long_cluster)
    amounts_per_cluster = []
    for cluster in longest_clusters:
        amounts = [0, 0]
        for read in cluster:
            if int(read) in sv_reads:
                amounts[0] += 1
            else:
                amounts[1] += 1
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
    print("Cluster purity score\t", purity_score)
    phasing_score = (phasing_sv_score + phasing_ref_score) / 2
    print("Phasing score\t\t", phasing_score)

    if len(clustering_matrix) == 2 and bad != True:
        print("clustering helemaal gelukt", clustering_matrix)

    purity_chart.append(int(purity_score))
    phasing_chart.append(int(phasing_score))
    return [int(purity_score), int(phasing_score)]


def make_clustering_matrix(matrix):
    clustering_matrix = {}
    for segment in range(len(matrix)):
        clustering_matrix[str(segment)] = {}
    for i in range(len(clustering_matrix)):
        for j in range(i + 1):
            mutations_in_common = 0
            if j == i:
                clustering_matrix[str(i)][str(j)] = 0
            else:
                amount_positions = len(matrix[i])
                for pos in range(len(matrix[i])):
                    if matrix[i][pos] is None or matrix[j][pos] is None:
                        amount_positions -= 1
                        continue
                    if matrix[i][pos] == -1 and matrix[j][pos] == -1:
                        amount_positions -= 1
                        continue
                    if matrix[i][pos] == matrix[j][pos]:
                        mutations_in_common += 1
                if amount_positions == 0:
                    clustering_matrix[str(i)][str(j)] = 0
                else:
                    clustering_matrix[str(i)][str(j)] = mutations_in_common / amount_positions
    return clustering_matrix