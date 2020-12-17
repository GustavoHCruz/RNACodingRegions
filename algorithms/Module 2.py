# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle   # Used for file operations
import re       # Used for regular expressions
# ======================================================

def verifySeq(PositionsSeq, ActualSeq, seqType):
    for i in PositionsSeq:
        if(i == ActualSeq):
            return seqType
    
    return "Neither"

# Open the result data from the last module
file = open("./assets/mod1.txt", "rb")
data = pickle.load(file) # The file that contains all the sequences generated in the module 1
# Declarations of variables that will be used in the algorithm
flag = True # Flag made to not allow empty sequences be included in the final result
combinations = [] # List that contains the final result, will append the other two lists
intronComb = []
exonComb = []
auxGT = 0 # Tells the position of the GT that we are looking for
auxAG = 0 # Tells the position of the AG that we are looking for
auxExternComb = [] # Contains all the possibilites for an entire sequence

for i in data: # For each sequence in data do
    auxExternComb = []
    auxGT = i[0].find("GT") # Find the first GT
    while(auxGT != -1): # While GT exists in sequence
        auxAG = i[0].find("AG",auxGT+2) # Find the first AG
        flag = False
        while(auxAG != -1): # While AG exists in sequence
            actualSeq = i[0][auxGT:auxAG+2]
            actualPositionSeq = [auxGT,auxAG+1]
            auxExternComb.append([actualSeq,actualPositionSeq,verifySeq(i[3],actualPositionSeq,"Intron")])
            auxAG = i[0].find("AG",auxAG+1) # Find the next AG in sequence, searching from the last AG position plus one
            flag = True

        auxGT = i[0].find("GT",auxGT+1) # Find the next GT in sequence, searching from the last GT position plus one

    intronComb.append(auxExternComb)

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
    
combinations = [exonComb,intronComb]

# Save the result data on mod2.txt archive
file = open("./assets/mod2.txt","wb")
pickle.dump(combinations,file)