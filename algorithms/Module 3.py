# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle               # Used for file operations
import sklearn_crfsuite     # Contains the functions to train, test and make predictions with the model
import random               # Used to generate random numbers
# ======================================================

# File generated in the previous module (you can edit the file name)
file = open("./assets/colletotrichumAtt_mod2.txt", "rb")
data = pickle.load(file)

# Samples is the variable that will contain all the sequences to be analyzed, according to the CRF input specification ([{'sequence': '...'}])
samples = []

# Labels is the variable that will contain all the sequences answers, for example, for samples[0], labels[0] can be or 'Exon', or 'Intron' or 'Neither'
labels = []

# The aux is the variable that will contain all the sequences of 'Neither'
aux = []

# The data[0] are all Exon sequences obtained in the last module. It is stored in the exons variable
exons = data[0]

# The data[1] are all Intron sequences obtained in the last module. It is stored in the introns variable
introns = data[1]

# Exon counter found
exonsCount = 0

# Intron counter found
intronsCount = 0

# While controller
i=0
while(i < len(exons)):                              # len(exons) = len(introns) = number of sequences obtained previously
    for seq in exons[i]:                            # for each subsequence in exons
        if(seq[-1] != "Neither"):                   # if the subsequence is not a 'Neither' one
            samples.append([{"sequence":seq[0]}])   # the sequence is appended to samples
            labels.append([seq[-1]])                # the sequence is appended to labels
            exonsCount += 1                         # exon counter is incremented
        else:
            aux.append([seq[0]])                    # The subsequence is a 'Neither' one and it's appended to aux

    for seq in introns[i]:                          # for each subsequence in introns
        if(seq[-1] != "Neither"):                   # if the subsequence is not a 'Neither' one
            samples.append([{"sequence":seq[0]}])   # the sequence is appended to samples
            labels.append([seq[-1]])                # the sequence is appended to labels
            intronsCount += 1                       # intron counter is incremented
        else:
            aux.append([seq[0]])                    # The subsequence is a 'Neither' one and it's appended to aux
    i+=1

# To balance the amount of 'Neither' in the base, the average between exons and introns is calculated
average = (exonsCount + intronsCount)/2

# For the average number of times, take a random aux element and append it to samples, putting 'Neither' in labels
for i in range(0,int(average)):
    randomNumber = random.randint(0,len(aux)-1)
    samples.append([{"sequence":aux[randomNumber][0]}])
    labels.append(["Neither"])
    aux.pop(randomNumber)

# Instance the CRF algorithm
clf = sklearn_crfsuite.CRF()

# Train the model
clf = clf.fit(samples, labels)

# Saves the trained model (you can edit the file name)
pickle.dump(clf, open("./results/colletotrichumAtt_model.sav", 'wb'))

# Saves the samples e labels
file = open("./assets/colletotrichumAtt_mod3.txt","wb")
pickle.dump([samples,labels],file)
# ========================================================================================================================