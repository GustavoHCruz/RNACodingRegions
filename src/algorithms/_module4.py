# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle               # Used for file operations
from copy import deepcopy

# Metrics is used to generate the final statistics on the trained model
from sklearn_crfsuite import metrics

# This import allows us to have a report of cross-validation results
from sklearn.model_selection import cross_val_predict
# ======================================================

# Import warnings filter
from warnings import simplefilter
# Ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)


# Load the model to be tested (you can edit the file name)
file = open("../../results/colletotrichum_model.sav", "rb")
clf = pickle.load(file)

# Load the data from module 3. You can edit the file name. Use the data corresponding to the model
file = open("../assets/colletotrichum_mod3.txt", "rb")
data = pickle.load(file)
trainSamples = data[0]
trainLabels = data[1]
testSamples = data[2]
testLabels = data[3]


pred = clf.predict(testSamples)
x = ['III', 'IIE', 'IEE', 'EEE', 'EEI', 'EII']
xFinal = ['Intron', 'Exon', 'Neither']
sorted_labels = sorted(x, key=lambda name: (name[1:], name[0]))

allSeqs = []

emptySeq = []
seq = []
currentSeq = 0

for i in range(len(testSamples)):
    if testSamples[i][0]['seqIndex'] == currentSeq:
        seq.append(pred[i][0])
    else:
        allSeqs.append(seq)
        seq = deepcopy(emptySeq)
        seq.append(pred[i][0])
        currentSeq +=1
allSeqs.append(seq)

finalPredictions = []
currentSeq = -1
for seq in allSeqs:
    currentSeq += 1

    seqLabel = ''
    seqLabelLength = len(seq)+2
    lastClass = 'I'
    lastClassStart = 0
    finalPrediction = {'id': currentSeq, 'seqLen': seqLabelLength, 'exons': [], 'introns': []}

    for i in range(seqLabelLength): # Surfacing through the seqLabel indexes
        if i == 0 or i == 1:
            label = seq[0][i]
            seqLabel += label
        if i == seqLabelLength-1 or i == seqLabelLength-2:
            label = seq[-1][i-(seqLabelLength-3)]
            seqLabel += label
        else:
            count = {'I': 0, 'E': 0}
            count[seq[i-2][2]] += 1
            count[seq[i-1][1]] += 1
            count[seq[i][0]] += 1
            label = 'I' if count['I'] > count['E'] else 'E'
            seqLabel += label
        
        if label != lastClass:
            if i > 0:
                if lastClass == 'I':
                    finalPrediction['introns'].append([lastClassStart, i-1])
                else:
                    finalPrediction['exons'].append([lastClassStart, i-1])
                lastClassStart = i
            lastClass = label

    if lastClass == 'I':
        finalPrediction['introns'].append([lastClassStart, seqLabelLength-1])
    else:
        finalPrediction['exons'].append([lastClassStart, seqLabelLength-1])

    finalPrediction['labelSeq'] = seqLabel
    finalPredictions.append(finalPrediction)

print(finalPredictions[37])


# print("Classification Report:")
# print(metrics.flat_classification_report(testLabels, pred, labels=sorted_labels, digits=2))
# ========================================================================================================================