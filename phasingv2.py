import os
import sys
# import plotly as py
# import plotly.graph_objs as go

matrix = []
posities = []
nucleotide_indexes = []
seqs = []


def preparation():
    with os.popen('/data/sambamba_v0.6.3  view ' + sys.argv[1]) as bam:
        read_number = 0
        sv_position = 167393746
        leftmost_mapping_position = None
        for line in bam:
            read = line.rstrip()
            columns = read.split("\t")
            cigar = columns[5]
            matrix.append([])
            nucleotide_indexes.append({})
            seqs.append(columns[9])
            cursor = int(columns[3]) - 1
            number = ""
            seq_base_index = 0
            if leftmost_mapping_position is None:
                leftmost_mapping_position = int(columns[3])
            elif leftmost_mapping_position > int(columns[3]):
                leftmost_mapping_position = int(columns[3])
            else:
                for none_position in range(leftmost_mapping_position, cursor + 1):
                    nucleotide_indexes[read_number][none_position] = None
            done = False
            for char in cigar:
                if done == True:
                    break
                if char.isdigit():
                    number += char
                elif char == "X":
                    for mismatch in range(int(number)):
                        cursor += 1
                        nucleotide_indexes[read_number][cursor] = seq_base_index
                        seq_base_index += 1
                    number = ""
                elif char == "=" or char == "D":
                    for i in range(int(number)):
                        cursor += 1
                        if char == "=":
                            nucleotide_indexes[read_number][cursor] = seq_base_index
                            seq_base_index += 1
                        elif char == "D":
                            nucleotide_indexes[read_number][cursor] = -1
                    number = ""
                elif char == "I":
                    seq_base_index += int(number)
                    number = ""
                elif char == "H":
                    number = ""
            if cursor < 167395746:
                for unknown_nuc_index in range(cursor, 167395746):
                    nucleotide_indexes[read_number][unknown_nuc_index] = None
            elif len(nucleotide_indexes[read_number]) < 1000:
                for unknown_nuc_index in range(sv_position, int(columns[3])):
                    nucleotide_indexes[read_number][unknown_nuc_index] = None
            read_number += 1


def create_na_values():
    for read in matrix:
        print(len(read))


def fill_matrix():
    with os.popen('/data/sambamba_v0.6.3  view ' + sys.argv[1]) as bam:
        read_number = 0
        teller = 0
        for line in bam:
            read = line.rstrip()
            columns = read.split("\t")
            cigar = columns[5]
            cursor = int(columns[3]) - 1
            number = ""
            done = False
            for char in cigar:
                if done == True:
                    break
                if char.isdigit():
                    number += char
                elif char == "X":
                    for mismatch in range(int(number)):
                        cursor += 1
                        if cursor >= 167395746:
                            done = True
                            break
                        if cursor >= 167393746:
                            if cursor not in posities:
                                teller += 1
                                idx = len(posities)
                                if len(posities) == 0:
                                    posities.append(cursor)
                                else:
                                    for index in range(len(posities)):
                                        if cursor < posities[index]:
                                            idx = posities.index(posities[index])
                                            posities.insert(idx, cursor)
                                            break
                                if len(posities) == idx:
                                    posities.append(cursor)
                                for read_number_bc in range(len(seqs)):
                                    if nucleotide_indexes[read_number_bc][cursor] == -1:
                                        matrix[read_number_bc].insert(idx, 0)
                                    elif nucleotide_indexes[read_number_bc][cursor] is None:
                                        matrix[read_number_bc].insert(idx, None)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'A':
                                        matrix[read_number_bc].insert(idx, 1)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'C':
                                        matrix[read_number_bc].insert(idx, 2)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'G':
                                        matrix[read_number_bc].insert(idx, 3)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'T':
                                        matrix[read_number_bc].insert(idx, 4)
                    number = ""
                elif char == "=" or char == "D":
                    for i in range(int(number)):
                        cursor += 1
                        if cursor >= 167395746:
                            done = True
                            break
                        if cursor >= 167393746:
                            if cursor not in posities and char == "D":
                                teller += 1
                                idx = len(posities)
                                if len(posities) == 0:
                                    posities.append(cursor)
                                else:
                                    for index in range(len(posities)):
                                        if cursor < posities[index]:
                                            idx = posities.index(posities[index])
                                            posities.insert(idx, cursor)
                                            break
                                if len(posities) == idx:
                                    posities.append(cursor)
                                for read_number_bc in range(len(seqs)):
                                    if nucleotide_indexes[read_number_bc][cursor] == -1:
                                        matrix[read_number_bc].insert(idx, 0)
                                    elif nucleotide_indexes[read_number_bc][cursor] is None:
                                        matrix[read_number_bc].insert(idx, None)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'A':
                                        matrix[read_number_bc].insert(idx, 1)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'C':
                                        matrix[read_number_bc].insert(idx, 2)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'G':
                                        matrix[read_number_bc].insert(idx, 3)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'T':
                                        matrix[read_number_bc].insert(idx, 4)
                    number = ""
                elif char == "I":
                    number = ""
                elif char == "H":
                    number = ""
            read_number += 1


def filter_matrix():
    to_be_deleted = []
    for pos in range(len(posities)):
        values_of_position = [0,0,0,0,0]
        for read in matrix:
            if read[pos] is not None:
                values_of_position[read[pos]] += 1
        true_calls = 0
        for value in values_of_position[1:]:
            if value > 0.25*len(matrix):
                true_calls += 1
        if true_calls < 2:
            to_be_deleted.append(posities[pos])
    for pos in to_be_deleted:
        index = posities.index(pos)
        for read in matrix:
            del read[index]
        del posities[index]


def heatmap():
    trace = go.Heatmap(z=matrix,
                       y=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                       x=posities)
    data = [trace]
    py.offline.plot(data, filename='labelled-heatmap')


def clustering():
    clustering_matrix = make_clustering_matrix()
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
        if highest_score < 0.3:
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
            merged_dict[str(j)] = sum_of_scores/2

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
    if len(clustering_matrix) > 2:
        print("Not able to phase to 2 clusters.")
        print("Started with:", len(seqs), "reads")
        print("Finished with:", len(clustering_matrix), "clusters")
    # print(len(clustering_matrix))
    # print(clustering_matrix)


def make_clustering_matrix():
    clustering_matrix = {}
    for read in range(len(matrix)):
        clustering_matrix[str(read)] = {}
    for i in range(len(clustering_matrix)):
        for j in range(i+1):
            mutations_in_common = 0
            if j == i:
                clustering_matrix[str(i)][str(j)] = 0
            else:
                amount_positions = len(matrix[i])
                for pos in range(len(matrix[i])):
                    if matrix[i][pos] is None or matrix[j][pos] is None:
                        amount_positions -= 1
                        continue
                    if matrix[i][pos] == matrix[j][pos]:
                        mutations_in_common += 1
                if amount_positions == 0:
                    clustering_matrix[str(i)][str(j)] = 0
                else:
                    clustering_matrix[str(i)][str(j)] = mutations_in_common / amount_positions
    return clustering_matrix


def main():
    preparation()
    fill_matrix()
    # create_na_values()
    filter_matrix()
    # heatmap()
    clustering()


main()

