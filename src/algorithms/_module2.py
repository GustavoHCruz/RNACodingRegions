# Functions ==============================================================================================================
# Integers List, List of Integers List and String -> String
# Check if some of the 'PositionsSeq' is equal to the 'CurrentSeq' and if true, 'seqType' is returned, otherwise, 'Neither' is returned
# Entry example: [84,93],[[0,61],[103,171],[227,318],[388,469]],'Intron'
# Output example: 'Neither'
def verifySeq(PositionsSeq, CurrentSeq, seqType):
    for i in PositionsSeq:
        if(i == CurrentSeq):
            return seqType
    return 'Neither'
assert(verifySeq([[0,61],[103,171],[227,318],[388,469]],[84,93],'Intron') == 'Neither')
assert(verifySeq([[0,22],[30,45],[80,100]],[30,45],'Intron') == 'Intron')
assert(verifySeq([[0,200]],[0,200],'Exon') == 'Exon')
assert(verifySeq([[0,200]],[5,30],'Exon') == 'Neither')


# 
def belongsToAExon(exons_list, index):
    for exon in exons_list:
        if index >= exon[0] and index <= exon[1]:
            return True
    return False
assert(belongsToAExon([[0,1],[2,3]], 0) == True)
assert(belongsToAExon([[0,1],[2,3]], 3) == True)
assert(belongsToAExon([[0,1],[2,3]], 4) == False)

# ========================================================================================================================

# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle   # Used for file operations
# ======================================================

# Open the result data from the last module (you can edit the file name)
file = open("../assets/colletotrichum_mod1.txt", "rb")

# The file that contains all the sequences generated in the module 1
data = pickle.load(file)

# Declarations of variables that will be used in the algorithm
flag = True         # Flag made to not allow empty sequences be included in the final result
combinations = []   # List that contains the final result, will append the other two lists
intronComb = []
exonComb = []
IIIComb = []
IIEComb = []
IEEComb = []
EEEComb = []
EEIComb = []
EIIComb = []
auxGT = 0           # Tells the position of the GT that we are looking for
auxAG = 0           # Tells the position of the AG that we are looking for
auxExternComb = []  # Contains all the possibilites for an entire sequence

# For each sequence in data do
for i in data:
    auxExternComb = []
    auxGT = i[0].find("GT")                     # Find the first GT
    while(auxGT != -1):                         # While GT exists in sequence
        auxAG = i[0].find("AG",auxGT+2)         # Find the first AG
        flag = False
        while(auxAG != -1):                     # While AG exists in sequence
            CurrentSeq = i[0][auxGT:auxAG+2]
            currentPositionSeq = [auxGT,auxAG+1]
            auxExternComb.append([CurrentSeq,currentPositionSeq,verifySeq(i[3],currentPositionSeq,"Intron")])
            auxAG = i[0].find("AG",auxAG+1)     # Find the next AG in sequence, searching from the last AG position plus one
            flag = True

        auxGT = i[0].find("GT",auxGT+1)         # Find the next GT in sequence, searching from the last GT position plus one

    intronComb.append(auxExternComb)




# Features:
#   Codon
#   Codon-1
#   Codon+1

# i = data[0] # Static
for i in data:
    seq = i[0]
    exons_list = i[1]

    # for j in range(len(seq)//3):
    #     codon = seq[3*j:3*(j+1)]
    #     codon_prev = seq[3*(j-1):3*j]
    #     codon_next = seq[3*(j+1):3*(j+2)]
    #     codon_label = ""
    #     for k in range(3*j,3*(j+1)):
    for j in range(len(seq)-2):
        codon = seq[j:j+3]
        codon_prev = seq[j-3:j]
        codon_next = seq[j+3:j+6]
        codon_label = ""
        for k in range(j,j+3):
            codon_label += "E" if belongsToAExon(exons_list, k) else "I"

        codonComb = [codon, codon_prev, codon_next, codon_label]

        if codon_label == "III":
            IIIComb.append(codonComb)
        elif codon_label == "IIE":
            IIEComb.append(codonComb)
        elif codon_label == "IEE":
            IEEComb.append(codonComb)
        elif codon_label == "EEE":
            EEEComb.append(codonComb)
        elif codon_label == "EEI":
            EEIComb.append(codonComb)
        elif codon_label == "EII":
            EIIComb.append(codonComb)

        # print(codon, codon_prev, codon_next, codon_label)


# For each sequence in data do
for i in data:
    auxExternComb = []    
    auxAG = 0
    flag = True
    while(flag):
        auxGT = i[0].find("GT",auxAG+2)
        while(auxGT != -1):
            CurrentSeq = i[0][auxAG:auxGT]
            currentPositionSeq = [auxAG,auxGT-1]
            auxExternComb.append([CurrentSeq,currentPositionSeq,verifySeq(i[1],currentPositionSeq,"Exon")])
            auxGT = i[0].find("GT",auxGT+1)
        CurrentSeq = i[0][auxAG:]
        currentPositionSeq = [auxAG,len(i[0])-1]
        auxExternComb.append([CurrentSeq,currentPositionSeq,verifySeq(i[1],currentPositionSeq,"Exon")])
        auxAG = i[0].find("AG",auxAG)
        if(auxAG == -1):
            flag = False
        else:
            auxAG += 2
    exonComb.append(auxExternComb)

# Append both exonComb and intronComb on combinations
combinations = [exonComb,intronComb,IIIComb,IIEComb,IEEComb,EEEComb,EEIComb,EIIComb]

# Save the result data (you can edit the file name)
file = open("../assets/colletotrichum_mod2.txt","wb")
pickle.dump(combinations,file)