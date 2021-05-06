import ProteinPrint

nucleotideCodeTable = {'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','C','G','T']}

def replaceLetters(sequence,result):
    counter = 0
    for letter in "RYSWKMBDHVN":
        if(sequence.find(letter) != -1):
            pos = 0
            for sub in nucleotideCodeTable[letter]:
                pos = sequence.find(letter, pos)
                replaceLetters(sequence[:pos] + sub + sequence[pos+1:],result)
        else:
            counter += 1
    if (counter == 11):
        result.append(sequence)

originalSeq = "TCCATCCCCAGACCGCAATCCGACCCCTTTAGGAGGTGATCGAGATCTGCGGTCGTCGGCCTTCGGAAAGGACATGTCGATTCTGACAACACCTCARCCTCCATTGTCGGTCGCCCTCGCCACCATGGGTATGTTTTTCTCCTCGCCCTCACTGCGGCATATTCCCCTATCCGCGCCGCGRCCGTCGCCTAACATGTGAATAGTATCATGATTGGTATGGGCC"

result = []

replaceLetters(originalSeq, result)

i = 1
for seq in result:
    print("\n\n============== " + str(i) + "º Sequence Being Analyzed ================")
    ProteinPrint.proteinTest(seq)
    i += 1