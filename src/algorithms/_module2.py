import random               # Used to generate random numbers

# Functions ==============================================================================================================

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

# Definitions
TEST_FRAC = 0.2

# Open the result data from the last module (you can edit the file name)
file = open("../assets/diaporthe_mod1.txt", "rb")

# The file that contains all the sequences generated in the module 1
data = pickle.load(file)

# Declarations of variables that will be used in the algorithm
combs = []
testCombs = []
testExonList = []
testIntronList = []

# Split the data between train and test
tenPercent = len(data) * TEST_FRAC
testData = []
for _ in range(int(tenPercent)):
    randomNumber = random.randint(0,len(data)-1)
    testData.append(data.pop(randomNumber))
    # testSamples.append(samples.pop(randomNumber))
    # testLabels.append(labels.pop(randomNumber))

# Features:
#   Codon
#   Codon-1
#   Codon+1

for i in range(len(data)):
    seq = data[i][0]
    exons_list = data[i][1]

    for j in range(len(seq)-2):
        codon = seq[j:j+3]
        codon_prev = seq[j-3:j]
        codon_next = seq[j+3:j+6]
        codon_label = ""
        for k in range(j,j+3):
            codon_label += "E" if belongsToAExon(exons_list, k) else "I"

        codonComb = [codon, codon_prev, codon_next, codon_label]

        combs.append(codonComb)

for i in range(len(testData)):
    seq = testData[i][0]
    exons_list = testData[i][1]
    introns_list = testData[i][3]

    testExonList.append(exons_list)
    testIntronList.append(introns_list)

    for j in range(len(seq)-2):
        codon = seq[j:j+3]
        codon_prev = seq[j-3:j]
        codon_next = seq[j+3:j+6]
        codon_label = ""
        for k in range(j,j+3):
            codon_label += "E" if belongsToAExon(exons_list, k) else "I"

        codonComb = [codon, codon_prev, codon_next, i, j, codon_label]

        testCombs.append(codonComb)


data = [combs, testCombs, testExonList, testIntronList]

# Save the result data (you can edit the file name)
file = open("../assets/diaporthe_mod2.txt","wb")
pickle.dump(data,file)