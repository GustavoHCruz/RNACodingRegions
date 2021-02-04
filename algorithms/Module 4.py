table = {'UUU':'F','CUU':'L','AUU':'I','GUU':'V','UUC':'F','CUC':'L','AUC':'I','GUC':'V','UUA':'L','CUA':'L','AUA':'I','GUA':'V','UUG':'L','CUG':'L','AUG':'M','GUG':'V','UCU':'S','CCU':'P','ACU':'T','GCU':'A','UCC':'S','CCC':'P','ACC':'T','GCC':'A','UCA':'S','CCA':'P','ACA':'T','GCA':'A','UCG':'S','CCG':'P','ACG':'T','GCG':'A','UAU':'Y','CAU':'H','AAU':'N','GAU':'D','UAC':'Y','CAC':'H','AAC':'N','GAC':'D','UAA':'Stop','CAA':'Q','AAA':'K','GAA':'E','UAG':'Stop','CAG':'Q','AAG':'K','GAG':'E','UGU':'C','CGU':'R','AGU':'S','GGU':'G','UGC':'C','CGC':'R','AGC':'S','GGC':'G','UGA':'Stop','CGA':'R','AGA':'R','GGA':'G','UGG':'W','CGG':'R','AGG':'R','GGG':'G'}

def translation(rna):
    protein = ""
    for start in range(0,len(rna),3):
        if table[rna[start:start+3]] != "Stop":
            protein += table[rna[start:start+3]]
    return protein

import pickle

file = open("./results/colletotrichum_model.sav", "rb")
clf = pickle.load(file)

seq = "CGGTGATGATGCGCCCAGAGCTGTCTTCCGTAAGTCTTCCCATCCGCAGACCGCAATCCGCCCCTTCAGGGGGGACTCAAATTTGCGGTCATCCAAACCGGGTGTGCTGTCGATACTAACCACCACGTAGCCTCCATTGTCGGTCGCCCTCGCCACCATGGGTATGTCTACTTCTCGCCCTCGCTGCGGTAATTTCCGCCCTCCGCGCCGCGATCTAACATGTGAATCAGTATCATGATTGGT"

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

for intron in intronComb:
    seq = seq.replace(intron[0],"")

seq = seq.replace("T","U")

while(True):
    if(len(seq)%3 != 0):
        seq = seq[0:-1]
    else:
        break

print(translation(seq))