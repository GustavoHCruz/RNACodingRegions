# Functions ==============================================================================================================
# Integers List, Lis of Integers List and String -> String
# Check if some of the 'PositionsSeq' is equal to the 'ActualSeq' and if true, 'seqType' is returned, otherwise, 'Neither' is returned
# Entry example: [84,93],[[0,61],[103,171],[227,318],[388,469]],'Intron'
# Output example: 'Neither'
def verifySeq(PositionsSeq, ActualSeq, seqType):
    for i in PositionsSeq:
        if(i == ActualSeq):
            return seqType
    return 'Neither'
assert(verifySeq([[0,61],[103,171],[227,318],[388,469]],[84,93],'Intron') == 'Neither')
assert(verifySeq([[0,22],[30,45],[80,100]],[30,45],'Intron') == 'Intron')
assert(verifySeq([[0,200]],[0,200],'Exon') == 'Exon')
assert(verifySeq([[0,200]],[5,30],'Exon') == 'Neither')
# ========================================================================================================================

# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle   # Used for file operations
# ======================================================

# Open the result data from the last module (you can edit the file name)
file = open("./assets/colletotrichumAtt_mod1.txt", "rb")

# The file that contains all the sequences generated in the module 1
data = pickle.load(file)

# Declarations of variables that will be used in the algorithm
flag = True         # Flag made to not allow empty sequences be included in the final result
combinations = []   # List that contains the final result, will append the other two lists
intronComb = []
exonComb = []
auxGT = 0           # Tells the position of the GT that we are looking for
auxAG = 0           # Tells the position of the AG that we are looking for
auxExternComb = []  # Contains all the possibilites for an entire sequence

# For each sequence in data do
for i in data:
    auxExternComb = []
    auxGT = i[0].find("GT")                     # Find the first GT
    while(auxGT != -1):                         # While GT exists in sequence
        auxAG = i[0].find("AG",auxGT+2)         # Find the first AG
        while(auxAG != -1):                     # While AG exists in sequence
            actualSeq = i[0][auxGT:auxAG+2]
            actualPositionSeq = [auxGT,auxAG+1]
            auxExternComb.append([actualSeq,actualPositionSeq,verifySeq(i[3],actualPositionSeq,"Intron")])
            auxAG = i[0].find("AG",auxAG+1)     # Find the next AG in sequence, searching from the last AG position plus one

        auxGT = i[0].find("GT",auxGT+1)         # Find the next GT in sequence, searching from the last GT position plus one

    intronComb.append(auxExternComb)

# For each sequence in data do
for i in data:
    auxExternComb = []    
    auxAG = 0
    flag = True
    while(flag):
        auxGT = i[0].find("GT",auxAG+2)
        while(auxGT != -1):
            actualSeq = i[0][auxAG:auxGT]
            actualPositionSeq = [auxAG,auxGT-1]
            auxExternComb.append([actualSeq,actualPositionSeq,verifySeq(i[1],actualPositionSeq,"Exon")])
            auxGT = i[0].find("GT",auxGT+1)
        actualSeq = i[0][auxAG:]
        actualPositionSeq = [auxAG,len(i[0])-1]
        auxExternComb.append([actualSeq,actualPositionSeq,verifySeq(i[1],actualPositionSeq,"Exon")])
        auxAG = i[0].find("AG",auxAG)
        if(auxAG == -1):
            flag = False
        else:
            auxAG += 2
    exonComb.append(auxExternComb)

# Append both exonComb and intronComb on combinations
combinations = [exonComb,intronComb]

# Save the result data (you can edit the file name)
file = open("./assets/colletotrichumAtt_mod2.txt","wb")
pickle.dump(combinations,file)