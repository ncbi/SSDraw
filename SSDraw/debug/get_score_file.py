import re

alignment = open("aligned.fasta").read().splitlines()

ID = "1oqy"
a_seq = ''
seq_found = 0
for i in alignment:
    if seq_found and i[0] == '>':
        break

    if i[0] == '>' and bool(re.search(ID.lower(), i.lower())):
        seq_found = 1
        continue

    if seq_found and i[0] != '>':
        a_seq += i


j = 0
for a in a_seq:
    if a != "-":
        k = j%11
        print(a+ " "+str(k))
        j+=1