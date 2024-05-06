import numpy as np
seq_name = []

def read_blosum(file_path):
    blosum_matrix = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        amino_acids = lines[0].split()
        for i in range(1, len(lines)):
            values = lines[i].split()
            for j in range(1, len(values)):
                blosum_matrix[(amino_acids[i - 1], amino_acids[j - 1])] = int(values[j])
    return blosum_matrix

def read_fasta(file_path):
    sequences = {}
    current_key = None
    current_sequence = ""

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if current_key is not None:
                    sequences[current_key] = list(current_sequence)

                current_key = line[1:]
                current_sequence = ""
                seq_name.append(current_key)
            else:
                current_sequence += line

        if current_key is not None:
            sequences[current_key] = list(current_sequence)

    return sequences

def global_alignment(seq1_key, seq2_key, sequences_dict, blosum_matrix, gap_penalty):
    seq1 = sequences_dict[seq1_key]
    seq2 = sequences_dict[seq2_key]

    m = len(seq1)
    n = len(seq2)

    dp_matrix = np.zeros((m + 1, n + 1))

    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(1, m + 1):
        dp_matrix[i, 0] = dp_matrix[i - 1, 0] + gap_penalty
        traceback_matrix[i, 0] = 1

    for j in range(1, n + 1):
        dp_matrix[0, j] = dp_matrix[0, j - 1] + gap_penalty
        traceback_matrix[0, j] = 2

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = dp_matrix[i - 1, j - 1] + blosum_matrix.get((seq1[i - 1], seq2[j - 1]))
            gap1_score = dp_matrix[i - 1, j] + gap_penalty
            gap2_score = dp_matrix[i, j - 1] + gap_penalty

            scores = [match_score, gap1_score, gap2_score]

            max_score_index = np.argmax(scores)
            dp_matrix[i, j] = scores[max_score_index]
            traceback_matrix[i, j] = max_score_index

    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n
    while i > 0 or j > 0:
        if traceback_matrix[i, j] == 0:  # Match or mismatch
            aligned_seq1.insert(0, seq1[i - 1])
            aligned_seq2.insert(0, seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 1:  # Gap in seq1
            aligned_seq1.insert(0, seq1[i - 1])
            aligned_seq2.insert(0, "-")
            i -= 1
        else:  # Gap in seq2
            aligned_seq1.insert(0, "-")
            aligned_seq2.insert(0, seq2[j - 1])
            j -= 1

    return aligned_seq1, aligned_seq2

def multiple_alignment(sequences,blosum,gap,tree_list):
    align_dict = {}
    seq1 = []
    seq2 = []
    alignseq1 = []
    alignseq2 = []
    temp1 = []
    temp2 = []
    count = 0
    for item in tree_list:

        items = item[0].split("-")
        if count == 0:
            align1,align2 = global_alignment(items[0],items[1],sequences,blosum,gap)
            for i in range(len(align2)):
                align1[i] = align1[i]+align2[i]
            str = items[0] + "," + items[1]
            align_dict[str] = align1
            count += 1

        # Multiple alignment
        else:
            m = 0
            n = 0

            if len(items[0].split(",")) == 1:
                seq1 = sequences[items[0]]
                m = len(seq1)
                temp1 = seq1
            elif len(items[0].split(",")) >= 1:
                alignseq1 = align_dict[items[0]]
                m = len(alignseq1)
                temp1 = alignseq1
            if len(items[1].split(",")) == 1:
                seq2 = sequences[items[1]]
                n = len(seq2)
                temp2 = seq2
            elif len(items[1].split(",")) >= 1:
                alignseq2 = align_dict[items[1]]
                n = len(alignseq2)
                temp2 = alignseq2

            dp_matrix = np.zeros((m + 1, n + 1))
            traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)
            for i in range(1, m + 1):
                dp_matrix[i, 0] = dp_matrix[i - 1, 0] + gap_penalty
                traceback_matrix[i, 0] = 1

            for j in range(1, n + 1):
                dp_matrix[0, j] = dp_matrix[0, j - 1] + gap_penalty
                traceback_matrix[0, j] = 2

            for i in range(1, m + 1):
                for j in range(1, n + 1):
                    match_score = 0
                    if len(seq1) > 0 and len(seq2)>0 and len(alignseq2) == 0 and len(alignseq1) == 0:
                        count = 0
                        break
                    elif len(seq1) > 0 and len(alignseq2)>0 and len(seq2) == 0 and len(alignseq1) == 0:
                        size = len(alignseq2[0])
                        match_score = dp_matrix[i - 1, j - 1]
                        for x in range(size):
                            if alignseq2[j-1][x] == "-":
                                match_score += gap_penalty
                            else:
                                match_score += blosum_matrix.get((seq1[i - 1], alignseq2[j-1][x]))


                    elif len(alignseq1) > 0 and len(seq2)>0 and len(seq1) == 0 and len(alignseq2) == 0:
                        size = len(alignseq1[0])
                        match_score = 0
                        for x in range(size):
                            if alignseq1[i - 1][x] == "-":
                                match_score += gap_penalty
                            else:
                                match_score += blosum_matrix.get((alignseq1[i - 1][x], seq2[j - 1]))
                    elif len(alignseq1) > 0 and len(alignseq2)>0 and len(seq2) == 0 and len(seq1) == 0:
                        size1 = len(alignseq1[0])
                        size2 = len(alignseq2[0])
                        match_score = 0
                        for x in range(size1):
                            for y in range(size2):
                                if alignseq1[i - 1][x] == "-" and alignseq2[j - 1][y] == "-":
                                    match_score += blosum_matrix.get(("*","*"))
                                elif (alignseq1[i - 1][x] == "-" and alignseq2[j - 1][y] != "-") or (alignseq1[i - 1][x] != "-" and alignseq2[j - 1][y] == "-"):
                                    match_score += gap_penalty
                                else:
                                    match_score += blosum_matrix.get((alignseq1[i - 1][x], alignseq2[j - 1][y]))


                    top_score = dp_matrix[i - 1, j] + gap_penalty
                    left_score = dp_matrix[i, j - 1] + gap_penalty
                    scores = [match_score, top_score, left_score]
                    max_score_index = np.argmax(scores)
                    dp_matrix[i, j] = scores[max_score_index]
                    traceback_matrix[i, j] = max_score_index
                    
                if count == 0:
                    break
            if count != 0:
                aligned_seq1 = []
                aligned_seq2 = []
                i, j = m, n

                while i > 0 or j > 0:
                    if traceback_matrix[i, j] == 0:  # Match or mismatch
                        aligned_seq1.insert(0, temp1[i - 1])
                        aligned_seq2.insert(0, temp2[j - 1])
                        i -= 1
                        j -= 1
                    elif traceback_matrix[i, j] == 1:  # Gap in seq1
                        char_ct = ""
                        aligned_seq1.insert(0, temp1[i - 1])
                        for x in range(len(temp2[0])):
                            char_ct += "-"
                        aligned_seq2.insert(0, char_ct)
                        i -= 1
                    else:  # Gap in seq2
                        char_ct = ""
                        for x in range(len(temp1[0])):
                            char_ct += "-"
                        aligned_seq1.insert(0, char_ct)
                        aligned_seq2.insert(0, temp2[j - 1])
                        j -= 1

                for i in range(len(aligned_seq2)):
                    aligned_seq1[i] = aligned_seq1[i] + aligned_seq2[i]
                str = items[0] + "," + items[1]
                align_dict[str] = aligned_seq1
                temp1.clear()
                temp2.clear()
            else:
                align1, align2 = global_alignment(items[0], items[1], sequences, blosum, gap)
                for i in range(len(align2)):
                    align1[i] = align1[i] + align2[i]
                str = items[0] + "," + items[1]
                align_dict[str] = align1
                count += 1
                
    str_list = str.split(",")
    for y in range(len(align_dict[str][0])):
        print("{:10}".format(str_list[y]), end="> ")
        for x in range(len(align_dict[str])):
            print(align_dict[str][x][y],end="")
        print()


def similarity_score(alignment):
    matched_count = 0
    aligned_length = 0

    min_length = min(len(alignment[0]), len(alignment[1]))

    for i in range(min_length):
        if alignment[0][i] == alignment[1][i] and alignment[0][i] != '-' and alignment[1][i] != '-':
            matched_count += 1
        if alignment[0][i] != '-' or alignment[1][i] != '-':
            aligned_length += 1

    if aligned_length == 0:
        return 0.0  # Alignment length is zero, similarity score should be zero

    similarity_score = matched_count / aligned_length
    return similarity_score


def create_similarity_matrix(sequences, blosum_matrix, gap_penalty):
    k = len(sequences)

    # Tüm olası çiftler arasında benzerlik skorlarını hesaplayın
    similarity_matrix = np.zeros((k, k))
    for i in range(k):
        similarity_matrix[i,i] = 1.0
    for i in range(k):
        for j in range(i + 1, k):
            alignment = global_alignment(seq_name[i], seq_name[j], sequences, blosum_matrix, gap_penalty)
            similarity_matrix[i, j] = similarity_score(alignment)
            similarity_matrix[j, i] = similarity_matrix[i, j]
    return similarity_matrix


# step2
def upgma(similarity_matrix, labels):
    n = len(labels)
    cluster_heights = {label: 0.0 for label in labels}

    cluster_list = []

    while n > 1:
        # Find the maximum similarity score in the matrix
        max_score = -1
        max_i, max_j = -1, -1
        for i in range(n):
            for j in range(i + 1, n):
                if similarity_matrix[i][j] > max_score:
                    max_score = similarity_matrix[i][j]
                    max_i, max_j = i, j

        # Calculate the new cluster label
        new_cluster_label = labels[max_i].replace("-", ",") + "-" + labels[max_j].replace("-", ",")

        # Update the cluster heights
        new_height = max_score / 2.0
        cluster_heights[new_cluster_label] = new_height

        # Update the similarity matrix and labels
        new_row = [(similarity_matrix[max_i][k] + similarity_matrix[max_j][k]) / 2.0 for k in range(n) if
                   k != max_i and k != max_j]
        new_col = [new_row[k] for k in range(n - 2)]

        similarity_matrix = [[similarity_matrix[i][j] for j in range(n) if j != max_i and j != max_j] for i in range(n)
                             if i != max_i and i != max_j]
        for i in range(n - 2):
            similarity_matrix[i].append(new_row[i])
        similarity_matrix.append(new_col + [0.0])

        labels = [labels[k] for k in range(n) if k != max_i and k != max_j]
        labels.append(new_cluster_label)  #rootu tutuyor

        cluster_list.append([new_cluster_label])

        # Reduce the number of clusters
        n -= 1

    return cluster_list


def find_max_element(matrix):
    max_value = 0
    max_i, max_j = -1, -1

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix[i])):
            if matrix[i][j] > max_value:
                max_value = matrix[i][j]
                max_i, max_j = i, j

    return max_i, max_j


fasta_file_path = "Input.txt"
blosum_file_path = "Blosum62.txt"
gap_penalty = int(input("Enter the gap penalty: "))
sequences = read_fasta(fasta_file_path)
blosum_matrix = read_blosum(blosum_file_path)


similarity_matrix = create_similarity_matrix(sequences, blosum_matrix, gap_penalty)

print("GUIDE TREE RESULT (STEP 2): ")
guide_tree = upgma(similarity_matrix.tolist(), seq_name)
print(guide_tree,"\n"), 

print("MULTIPLE ALIGNMENT RESULT (STEP 3): ")
multiple_alignment(sequences,blosum_matrix,gap_penalty,guide_tree)


