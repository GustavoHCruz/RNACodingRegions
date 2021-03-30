# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle               # Used for file operations
import sklearn_crfsuite     # Contains the functions to train, test and make predictions with the model
# ======================================================


# File generated in the previous module (you can edit the file name)
file = open("../assets/diaporthe_mod2.txt", "rb")
data = pickle.load(file)

# Samples is the variable that will contain all the sequences to be analyzed, according to the CRF input specification ([{'sequence': '...'}])
samples = []
testSamples = []

# Labels is the variable that will contain all the sequences answers, for example, for samples[0], labels[0] can be or 'Exon', or 'Intron' or 'Neither'
labels = []
testLabels = []


combs = data[0]
testCombs = data[1]
testExonList = data[2]
testIntronList = data[3]

# testCombs = sorted(testCombs, key=lambda k: (k[3], k[4]))

for codon in combs:
    samples.append([{"codon": codon[0], "codon-1": codon[1], "codon+1": codon[2]}])
    labels.append([codon[-1]])

for codon in testCombs:
    testSamples.append([{"codon": codon[0], "codon-1": codon[1], "codon+1": codon[2], "seqIndex": codon[3], "codonIndex": codon[4]}])
    testLabels.append([codon[-1]])



# Instance the CRF algorithm
clf = sklearn_crfsuite.CRF()

# Train the model
clf = clf.fit(samples, labels)

# Saves the trained model (you can edit the file name)
pickle.dump(clf, open("../../results/diaporthe_model.sav", 'wb'))

# Saves the samples e labels
file = open("../assets/diaporthe_mod3.txt","wb")
pickle.dump([testSamples, testLabels, testExonList, testIntronList], file)
# ========================================================================================================================