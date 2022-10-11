import numpy as np
from Bio import SeqIO

START_CODON = 'ATG'
STOP_CODONS = ['TAA', 'TAG', 'TGA']
LETTERS = 'ACGT'


def get_codon_frames(seq):
    frames = [[seq[i:i + 3] for i in range(0, len(seq), 3)], [seq[i:i + 3] for i in range(1, len(seq), 3)],
              [seq[i:i + 3] for i in range(2, len(seq), 3)],
              [seq.reverse_complement()[i:i + 3] for i in range(0, len(seq), 3)],
              [seq.reverse_complement()[i:i + 3] for i in range(1, len(seq), 3)],
              [seq.reverse_complement()[i:i + 3] for i in range(2, len(seq), 3)]]
    return frames


def find_start_stop_pairs(seq_record):
    i = 0
    start_stop_pair = []
    while i < len(seq_record):
        if seq_record[i] == START_CODON:
            start_pos = j = i
            while j < len(seq_record):
                if seq_record[j] in STOP_CODONS:
                    end_pos = i = j
                    start_stop_pair.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    break
                j += 1
        i += 1
    return start_stop_pair


def find_farthest_codon_for_each_stop_codon(seq_record):
    i = 0
    codon_list = []
    while i < len(seq_record):
        if seq_record[i] in STOP_CODONS:
            start_pos = i
            j = i + 1
            end_pos = -1
            while j < len(seq_record):
                if seq_record[j] == START_CODON:
                    end_pos = j
                elif seq_record[j] in STOP_CODONS or j + 1 == len(seq_record):
                    if end_pos != -1:
                        codon_list.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codon_list


def find_codon_frequency(codon_list):
    alphabet_list = []
    for c1 in LETTERS:
        for c2 in LETTERS:
            for c3 in LETTERS:
                codon = c1 + c2 + c3
                alphabet_list.append(codon)
    codon_alphabet = dict.fromkeys(alphabet_list, 0.0)

    total = 0
    for orf in codon_list:
        codons_in_orf = [orf[i:i + 3] for i in range(0, len(orf), 3)]
        total += len(codons_in_orf)
        for codon in codons_in_orf:
            codon_alphabet[codon] += 1
    for key in codon_alphabet:
        codon_alphabet[key] = (codon_alphabet[key] / total) * 100
    return codon_alphabet


def find_dicodon_frequency(dicodon_list):
    alphabet_list = []
    for c1 in LETTERS:
        for c2 in LETTERS:
            for c3 in LETTERS:
                for c4 in LETTERS:
                    for c5 in LETTERS:
                        for c6 in LETTERS:
                            dicodon = c1 + c2 + c3 + c4 + c5 + c6
                            alphabet_list.append(dicodon)
    dicodon_alphabet = dict.fromkeys(alphabet_list, 0.0)

    total = 0
    for orf in dicodon_list:
        dicodons_in_orf = [orf[i:i + 6] for i in range(0, len(orf), 6)]
        total += len(dicodons_in_orf)
        for dicodon in dicodons_in_orf:
            if len(dicodon) == 6:
                dicodon_alphabet[dicodon] += 1
    for key in dicodon_alphabet:
        dicodon_alphabet[key] = (dicodon_alphabet[key] / total) * 100
    return dicodon_alphabet


def print_distance_matrix(frequencies):
    names = [key for key in frequencies]
    alphabet = frequencies[names[0]].keys()

    print("*******", end=' ')
    for name in names:
        print(name, end=' ')

    for i in range(0, len(frequencies)):
        print('\n' + names[i], end=' ')
        for j in range(0, len(frequencies)):
            total_points = sum(list(map(lambda k: abs(frequencies[names[i]][k] - frequencies[names[j]][k]), alphabet)))
            print("%.2f" % total_points, end=' ')


if __name__ == '__main__':
    files = ["bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta", "mamalian1.fasta",
             "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"]

    codon_frequencies = {}
    dicodon_frequencies = {}

    for file in files:
        record = SeqIO.read("data/" + file, "fasta")
        codon_frames = get_codon_frames(record.seq)

        print("*****" + record.id + "*****")
        # 1.
        # print("1. All codons in a sequence:")
        start_stop_pairs = np.concatenate(list(map(find_start_stop_pairs, codon_frames)))
        # print(start_stop_pairs)

        # 2.
        # print("2. Each stop codons farthest start codon")
        farthest_start_codons = list(map(find_farthest_codon_for_each_stop_codon, codon_frames))
        # print(farthest_start_codons)

        # 3.
        longest_fragments = list(filter(lambda orf: len(orf) > 100, start_stop_pairs))
        # print("3. Fragments that have a length more than 100 symbols:")
        # print(longest_fragments)

        # 4.
        # print("4. Codon frequency (%):")
        codon_frequencies[record.id] = find_codon_frequency(longest_fragments)
        # print(find_codon_frequency(codons))

        # print("4. Dicodon frequency (%):")
        dicodon_frequencies[record.id] = find_dicodon_frequency(longest_fragments)
        # print(find_dicodon_frequency(codons))

    print("\nCodon frequency matrix:")
    print_distance_matrix(codon_frequencies)

    print("\n\nDicodon frequency matrix:")
    print_distance_matrix(dicodon_frequencies)
