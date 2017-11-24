import os
import sys
import plotly as py
import plotly.graph_objs as go
from sklearn.cluster import KMeans
import numpy as np

matrix = []
posities = []
nucleotide_indexes = []
seqs = []

def preparation():
    with os.popen('/data/sambamba_v0.6.3  view ' + sys.argv[1]) as bam:
        read_number = 0
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
            done = False
            for char in cigar:
                if done == True:
                    break
                if char.isdigit():
                    number += char
                elif char == "X":
                    for mismatch in range(int(number)):
                        cursor+=1
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
            read_number += 1

def main():
    preparation()
    with os.popen('/data/sambamba_v0.6.3  view ' + sys.argv[1]) as bam:
        read_number = 0
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
                        cursor+=1
                        if cursor >= 149006000:
                            # print(cursor, columns[9][seq_base_index])
                            if cursor not in posities:
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
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'A':
                                        matrix[read_number_bc].insert(idx, 1)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'C':
                                        matrix[read_number_bc].insert(idx, 2)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'G':
                                        matrix[read_number_bc].insert(idx, 3)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'T':
                                        matrix[read_number_bc].insert(idx, 4)
                        if cursor >= 149006955:
                            done = True
                            break
                    number = ""
                elif char == "=" or char == "D":
                    for i in range(int(number)):
                        cursor += 1
                        if cursor >= 149006000:
                            if cursor in posities and char == "D":
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
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'A':
                                        matrix[read_number_bc].insert(idx, 1)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'C':
                                        matrix[read_number_bc].insert(idx, 2)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'G':
                                        matrix[read_number_bc].insert(idx, 3)
                                    elif seqs[read_number_bc][nucleotide_indexes[read_number_bc][cursor]] == 'T':
                                        matrix[read_number_bc].insert(idx, 4)
                        if cursor >= 149006955:
                            done = True
                            break
                    number = ""
                elif char == "I":
                    number = ""
                elif char == "H":
                    number = ""
            read_number += 1
#
main()
# for rd in matrix:
#     print(len(rd))
# heatmap_pos = []
# for coordinate in range(0,956):
#     heatmap_pos.append(coordinate)

# print("\n\n", nucleotide_indexes[13][149006000], matrix[13][0])

# trace = go.Heatmap(z = matrix,
#                    y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
#                    x = posities)
# data = [trace]
# py.offline.plot(data, filename = 'labelled-heatmap')

#,
#                    colorscale=[[0, 'rgb(128, 128, 128)'], [1, 'rgb(0, 179, 0)'], [2, 'rgb(0, 0, 255)'],
#                                [3, 'rgb(255, 204, 0)'], [4, 'rgb(255, 0, 0)']]


X = np.array(matrix)
kmeans = KMeans(n_clusters=2, random_state=0)
print(kmeans.fit_predict(X))
