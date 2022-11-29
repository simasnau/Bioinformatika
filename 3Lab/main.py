from bioinfokit import analys
import matplotlib.pyplot as plt
from collections import Counter

analys.format.fq_qual_var(file='reads_for_analysis.fastq')


read_file = open('reads_for_analysis.fastq', 'r')
lines = read_file.readlines()
read_file.close()

sequences = {}

i = 0

lines_count = len(lines)
while i < lines_count:
    name = lines[i]
    seq = lines[i+1]
    sequences[name] = seq
    i=i+4
    if i >= lines_count:
        break

percentages = {}

for name in sequences:
    seq = sequences[name]

    seq_length = len(seq)
    occurences = 0
    for c in seq:
        if c == 'C' or c == 'G': 
            occurences = occurences + 1
    percentage = int((occurences / seq_length) * 100)

    percentages[name] = percentage

percentages_list = list(percentages.values())

binwidth = 1

plt.hist(percentages_list, bins=range(min(percentages_list), max(percentages_list) + binwidth, binwidth))
plt.show()

peaks = [34, 50, 69]

for peak in peaks:
    listOfKeys = [key  for (key, value) in percentages.items() if value == peak]
    print('5 Sequences for', peak)

    keys = listOfKeys[:5]
    for key in keys:
        print(key,':', sequences[key])
