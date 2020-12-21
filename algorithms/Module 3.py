import pickle
import sklearn_crfsuite
import random

# from sklearn_crfsuite import metrics
# from sklearn.model_selection import cross_val_predict

file = open("./assets/mod2.txt", "rb")
data = pickle.load(file)

samples = []
labels = []
aux = []

exons = data[0]
introns = data[1]
exonsCount = 0
intronsCount = 0

i=0
while(i < len(exons)):
    for seq in introns[i]:
        if(seq[-1] != "Neither"):
            samples.append([{"sequence":seq[0]}])
            intronsCount += 1
            labels.append([seq[-1]])
        
        aux.append([seq[0]])
    for seq in exons[i]:
        if(seq[-1] != "Neither"):
            samples.append([{"sequence":seq[0]}])
            exonsCount += 1
            labels.append([seq[-1]])

        aux.append([seq[0]])
    i+=1

average = (exonsCount + intronsCount)/2
for i in range(0,int(average)):
    randomNumber = random.randint(0,len(aux)-1)
    samples.append([{"sequence":aux[randomNumber]}])
    labels.append(["Neither"])
    aux.pop(randomNumber)

clf = sklearn_crfsuite.CRF(algorithm="lbfgs", c1=0.07, c2=0.09, all_possible_transitions=True)
clf = clf.fit(samples, labels)

# pred = cross_val_predict(clf, samples, labels, cv=10)
# print("Classification_report:")
# labels = ['Intron', 'Exon', 'Neither']
# metrics.flat_f1_score(labels, pred, average='weighted', labels=labels)
# sorted_labels = sorted(labels, key=lambda name: (name[1:], name[0]))
# print(metrics.flat_classification_report(
#     labels, pred, labels=sorted_labels, digits=2))

pickle.dump(clf, open("./results/model.sav", 'wb'))