import ProteinPrint

nucleotideCodeTable = {'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','C','G','T']}

def replaceLetters(sequence, initialPosition):
    appendFlag = True
    for i in range(initialPosition,len(sequence)):
        letter = sequence[i]
        if (letter in "RYSWKMBDHVN"):
            appendFlag = False
            for possibility in nucleotideCodeTable[letter]:
                newSeq = sequence[:i] + possibility + sequence[i+1:]
                replaceLetters(newSeq, i+1)
            break

    if(appendFlag):
        result.append(sequence)
        
originalSeq = "TCCATTCCCAGACCGCAATCCCACCCCTTTAGGAGGTGATCRAGATCTGCSGTCGKCGGCTTTCASAAAGGACATGTCRATTCTGACAATACCTCAGCCTCCATTGTCGGTCGCCCTCGCCACCATGGGTATGTTTTTCTCCTCKCCCTCACTGCGGSATATTCCCCTATTCGCGCCGCGACCGTCSCCTAACATGTGAATAGTATCATGATTGGTATGGGCC"

result = []

replaceLetters(originalSeq, 0)

i = 1
for seq in result:
    print("\n\n============== " + str(i) + "º Sequence Being Analyzed ================")
    ProteinPrint.proteinTest(seq)
    i += 1