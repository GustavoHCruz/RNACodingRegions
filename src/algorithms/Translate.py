import re

f = open("./assets/coddon_table.txt")

table = {}

for line in f:
    line = line.strip()
    line = [ (codon.split(" ")[0], codon.split(" ")[1]) for codon in re.split(r"[ ]{2,3}", line)if codon]
    for r in line:
        table[r[0]] = r[1]

print(table)