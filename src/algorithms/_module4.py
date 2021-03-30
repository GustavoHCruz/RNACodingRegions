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

def findExonsAndIntrons(seq: str):
    lastClass = 'I'
    lastClassStart = 0
    exons = []
    introns = []
    for i in range(len(seq)):
        label = seq[i]
        if label != lastClass:
            if i > 0:
                if lastClass == 'I':
                    introns.append([lastClassStart, i-1])
                else:
                    exons.append([lastClassStart, i-1])
                lastClassStart = i
            lastClass = label
    if lastClass == 'I':
        introns.append([lastClassStart, len(seq)-1])
    else:
        exons.append([lastClassStart, len(seq)-1])
    return exons, introns


def aliasing(seq: str, rule: int):
    newSeq = seq[:]
    noises = {}
    for i in range(len(newSeq)+1):
        if i == 0:
            counter = 1
            lastClass = newSeq[0]
            lastClassStart = 0
            prevClassStart = -1
        elif i == len(newSeq):
            if counter <= rule:
                dictAux = {}
                dictAux['counter'] = counter
                dictAux['start'] = lastClassStart
                dictAux['end'] = i
                dictAux['class'] = lastClass
                dictAux['prev'] = prevClassStart
                noises[lastClassStart] = dictAux
        else:
            if newSeq[i] == lastClass:
                counter += 1
            else:
                if counter <= rule:
                    dictAux = {}
                    dictAux['counter'] = counter
                    dictAux['start'] = lastClassStart
                    dictAux['end'] = i
                    dictAux['class'] = lastClass
                    dictAux['prev'] = prevClassStart
                    noises[lastClassStart] = dictAux
                prevClassStart = lastClassStart
                lastClass = newSeq[i]
                lastClassStart = i
                counter = 1
    # print(noises.keys())
    # print(newSeq)
    while noises:
        reg = noises.pop(min(noises.items(), key=lambda item: item[1]['counter'])[0])
        # print(noises.keys())
        # print(reg)
        newSeq = newSeq[:reg['start']] + reg['counter'] * ('E' if reg['class'] == 'I' else 'I') + newSeq[reg['end']:]
        # print(newSeq, "oi")
        prev = reg['prev']
        if prev in noises.keys():
            noises[prev]['counter'] += reg['counter']
        cond = True
        current = reg
        while cond:
            cond = False
            for noise in noises.items():
                if noise[1]['prev'] == current['start']:
                    if noise[1]['class'] == reg['class']:
                        newSeq = newSeq[:noise[1]['start']] + noise[1]['counter'] * ('E' if noise[1]['class'] == 'I' else 'I') + newSeq[noise[1]['end']:]
                    if prev in noises.keys():
                        noises[prev]['counter'] += noise[1]['counter']
                    current = noise[1]
                    noises.pop(noise[0])
                    # print(noise[1])
                    # print(noises.keys())
                    # print(newSeq, "tchau")
                    cond = True
                    break
        if prev in noises.keys() and noises[prev]['counter'] > rule:
            noises.pop(prev)
    # print(newSeq)
    return newSeq

assert(aliasing('IIIEIEIEEEEEEEIIEEEEEE', 5) == 'IIIIIIIEEEEEEEEEEEEEEE')
assert(aliasing('EEEEEEEEIIIIIIIIEIEEEEEEIIIEIEE', 5) == 'EEEEEEEEIIIIIIIIIIEEEEEEIIIIIII')
assert(aliasing('EEEEEEEEIIIIIIIIEIEEEEEEIIIEIEE', 4) == 'EEEEEEEEIIIIIIIIIIEEEEEEIIIIIII')


def indexToLabel(exons: list, index: int):
    for exon in exons:
        if index >= exon[0] and index <= exon[1]:
            return 'Exon'
    return 'Intron'


def findSplicingRegs(exonList: list, intronList: list):
    auxList = deepcopy(exonList)
    auxList.extend(intronList)
    auxDict = {}
    splicingRegs = []
    for x in auxList:
        auxDict[x[0]] = x[1]
    i = 0
    while auxDict[i]+1 in auxDict.keys():
        i = auxDict[i]+1
        splicingRegs.append(i-1)
    return splicingRegs

assert(findSplicingRegs([[0,1], [7,9]], [[2,6], [10,15]]) == [1,6,9])


def listsCompare(list1: list, list2: list):
    matches = [x for x in list1 if x in list2]
    misses = [x for x in list1 if x not in list2]
    skips = [x for x in list2 if x not in list1]
    distance = max(len(misses), len(skips))
    return len(matches), distance

assert(listsCompare([1,6,9], [1,8,9]) == (2, 1))


def calcMetrics(finalPredictions: list, testExonList: list, testIntronList: list):
    total = len(finalPredictions)
    realLabels = []
    predLabels = []
    yFinal = ['Intron', 'Exon']
    sorted_labels = sorted(yFinal)
    matchesSum = 0
    distanceSum = 0
    for i in range(total):
        for j in range(len(finalPredictions[i]['labelSeq'])):
            predLabels.append(['Intron' if finalPredictions[i]['labelSeq'][j] == 'I' else 'Exon'])
            realLabels.append([indexToLabel(testExonList[i], j)])
        matches, distance = listsCompare(findSplicingRegs(testExonList[i], testIntronList[i]), findSplicingRegs(finalPredictions[i]['exons'], finalPredictions[i]['introns']))
        matchesSum += matches
        distanceSum += distance
    splicingPrecision = matchesSum / (matchesSum + distanceSum)    

    print("Final Classification Report:")
    print(metrics.flat_classification_report(realLabels, predLabels, labels=sorted_labels, digits=2))
    print("\nSplicing Precision: {:.3f}" .format(splicingPrecision))




# Load the model to be tested (you can edit the file name)
file = open("../../results/colletotrichum_model.sav", "rb")
clf = pickle.load(file)

# Load the data from module 3. You can edit the file name. Use the data corresponding to the model
file = open("../assets/colletotrichum_mod3.txt", "rb")
data = pickle.load(file)
testSamples = data[0]
testLabels = data[1]
testExonList = data[2]
testIntronList = data[3]

# print(len(testSamples))
# print(len(testLabels))
# print(len(testExonList))
# print(len(testIntronList))


pred = clf.predict(testSamples)

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
    finalPrediction = {'id': currentSeq, 'seqLen': seqLabelLength}

    for i in range(seqLabelLength): # Surfacing through the seqLabel indexes
        if i == 0 or i == 1:
            label = seq[0][i]
            seqLabel += label
        elif i == seqLabelLength-1 or i == seqLabelLength-2:
            label = seq[-1][i-(seqLabelLength-3)]
            seqLabel += label
        else:
            count = {'I': 0, 'E': 0}
            count[seq[i-2][2]] += 1
            count[seq[i-1][1]] += 1
            count[seq[i][0]] += 1
            label = 'I' if count['I'] > count['E'] else 'E'
            seqLabel += label
        
    finalPrediction['labelSeq'] = aliasing(seqLabel, 6)
    finalPrediction['exons'], finalPrediction['introns'] = findExonsAndIntrons(finalPrediction['labelSeq'])
    finalPredictions.append(finalPrediction)

calcMetrics(finalPredictions, testExonList, testIntronList)
# print(finalPredictions[10])
# finalPredictions[10]['labelSeq'] = aliasing(finalPredictions[10]['labelSeq'], 4)
# print(finalPredictions[10])
# print(findExonsAndIntrons(finalPredictions[10]['labelSeq']))
# print(testExonList[10], testIntronList[10])
# print(testIntronList[10])


# ========================================================================================================================