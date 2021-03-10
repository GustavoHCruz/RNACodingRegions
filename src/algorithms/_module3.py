# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle               # Used for file operations
import sklearn_crfsuite     # Contains the functions to train, test and make predictions with the model
# ======================================================


# File generated in the previous module (you can edit the file name)
file = open("../assets/colletotrichum_mod2.txt", "rb")
data = pickle.load(file)

# Samples is the variable that will contain all the sequences to be analyzed, according to the CRF input specification ([{'sequence': '...'}])
samples = []
testSamples = []

# Labels is the variable that will contain all the sequences answers, for example, for samples[0], labels[0] can be or 'Exon', or 'Intron' or 'Neither'
labels = []
testLabels = []


combs = data[0]
testCombs = data[1]

testCombs = sorted(testCombs, key=lambda k: (k[3], k[4]))

for codon in combs:
    samples.append([{"codon": codon[0], "codon-1": codon[1], "codon+1": codon[2]}])
    labels.append([codon[-1]])

for codon in testCombs:
    testSamples.append([{"codon": codon[0], "codon-1": codon[1], "codon+1": codon[2], "seqIndex": codon[3], "codonIndex": codon[4]}])
    testLabels.append([codon[-1]])


# # Exon counter found
# exonsCount = 0

# # Intron counter found
# intronsCount = 0


# # While controller
# i=0
# while(i < len(exons)):                              # len(exons) = len(introns) = number of sequences obtained previously
#     for seq in exons[i]:                            # for each subsequence in exons
#         if(seq[-1] != "Neither"):                   # if the subsequence is not a 'Neither' one
#             samples.append([{"sequence":seq[0]}])   # the sequence is appended to samples
#             labels.append([seq[-1]])                # the sequence is appended to labels
#             exonsCount += 1                         # exon counter is incremented
#         else:
#             aux.append([seq[0]])                    # The subsequence is a 'Neither' one and it's appended to aux

#     for seq in introns[i]:                          # for each subsequence in introns
#         if(seq[-1] != "Neither"):                   # if the subsequence is not a 'Neither' one
#             samples.append([{"sequence":seq[0]}])   # the sequence is appended to samples
#             labels.append([seq[-1]])                # the sequence is appended to labels
#             intronsCount += 1                       # intron counter is incremented
#         else:
#             aux.append([seq[0]])                    # The subsequence is a 'Neither' one and it's appended to aux
#     i+=1


    

# print(samples)
# print(labels)
# print()
# print(testSamples)
# print(testLabels)

# # To balance the amount of 'Neither' in the base, the average between exons and introns is calculated
# average_train = (exonsCount + intronsCount)/2
# average_test = average_train * TEST_FRAC

# # For the average number of times, take a random aux element and append it to samples, putting 'Neither' in labels
# for i in range(int(average_train)):
#     randomNumber = random.randint(0,len(aux)-1)
#     samples.append([{"sequence":aux[randomNumber]}])
#     labels.append(["Neither"])
#     aux.pop(randomNumber)

# for i in range(int(average_test)):
#     randomNumber = random.randint(0,len(aux)-1)
#     testSamples.append([{"sequence":aux[randomNumber]}])
#     testLabels.append(["Neither"])
#     aux.pop(randomNumber)

# Instance the CRF algorithm
clf = sklearn_crfsuite.CRF()

# Train the model
clf = clf.fit(samples, labels)

# Saves the trained model (you can edit the file name)
pickle.dump(clf, open("../../results/colletotrichum_model.sav", 'wb'))

# Saves the samples e labels
file = open("../assets/colletotrichum_mod3.txt","wb")
pickle.dump([samples,labels,testSamples,testLabels],file)
# ========================================================================================================================