table = {'UUU':'F','CUU':'L','AUU':'I','GUU':'V','UUC':'F','CUC':'L','AUC':'I','GUC':'V','UUA':'L','CUA':'L','AUA':'I','GUA':'V','UUG':'L','CUG':'L','AUG':'M','GUG':'V','UCU':'S','CCU':'P','ACU':'T','GCU':'A','UCC':'S','CCC':'P','ACC':'T','GCC':'A','UCA':'S','CCA':'P','ACA':'T','GCA':'A','UCG':'S','CCG':'P','ACG':'T','GCG':'A','UAU':'Y','CAU':'H','AAU':'N','GAU':'D','UAC':'Y','CAC':'H','AAC':'N','GAC':'D','UAA':'Stop','CAA':'Q','AAA':'K','GAA':'E','UAG':'Stop','CAG':'Q','AAG':'K','GAG':'E','UGU':'C','CGU':'R','AGU':'S','GGU':'G','UGC':'C','CGC':'R','AGC':'S','GGC':'G','UGA':'Stop','CGA':'R','AGA':'R','GGA':'G','UGG':'W','CGG':'R','AGG':'R','GGG':'G'}

def translation(rna):
    protein = ""
    for start in range(0,len(rna),3):
        if (len(rna[start:start+3]) == 3 and table[rna[start:start+3]] != "Stop"):
            protein += table[rna[start:start+3]]
    return protein

import pickle

file = open("./results/diaporthe_model.sav", "rb")
clf = pickle.load(file)

originalSeq = "CGGTGCTGCTTTCTGGTGCGTACCASCTCCAGCTCCGAGCCTACCACCGCGATGATCGACGCGCGACAAGGCGAGCTCGAAGCATCGATACTGACCTCGGTTCTTTAGGCAAACAATCTCTGGCGAGCACGGTCTCAACAGCAATGGSRTGTATGTACCTCCTATTCCCTGACTACTGACCTCGGCCTCTCCTCCGGSTTGGGACTGACGATSGCACAGTTACAACGGCACTTCCGAACTCCAACTCGAGCGCATGAACATCTACTTCAACGAGGTAAGTCAATAGCCACGTCGTCRATTCKAGTTTGACCCTCTCGGCATGGGTGACTGCCGCCGCCAAACCCTTGCTAACGCGTTCTCGCCCAGGCCTCCGGCAACAAGTATGTGCCCCGCGCCGTCCTCGTCGATCTCGA"


originalSeq = originalSeq.upper()
seq = originalSeq

flag = True
intronComb = []
auxGT = 0
auxAG = 0

intronComb = []
auxGT = seq.find("GT")
while(auxGT != -1):
    auxAG = seq.find("AG",auxGT+2)
    flag = False
    while(auxAG != -1):
        actualSeq = seq[auxGT:auxAG+2]
        if(['Intron'] == clf.predict_single([{'sequence':actualSeq}])):
            intronComb.append([actualSeq])

        auxAG = seq.find("AG",auxAG+1)
        flag = True

    auxGT = seq.find("GT",auxGT+1)

print("Foram encontrados",len(intronComb),"Introns na sequência:")

for intron in intronComb:
    pos = seq.find(intron[0])
    print("Posição do Intron na sequência:",pos,"-",pos+len(intron[0]))
    print("Intron:",intron[0])
for intron in intronComb:
    seq = seq.replace(intron[0],"")

aux = originalSeq
exonComb = []
for intron in intronComb:
    aux = aux.split(intron[0])
    exonComb.append(aux[0])
    aux = aux[1]
exonComb.append(aux)
print("\nForam encontrados",len(exonComb),"Exons na sequência:")
for exon in exonComb:
    print("Exon:",exon)

print("\nSequência pós splicing:",seq)
seq = seq.replace("T","U")
print("\nTradução:",translation(seq))
print("Tradução +1:",translation(seq[1:]))
print("Tradução +2:",translation(seq[2:]))