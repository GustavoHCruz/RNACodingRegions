import pickle
import sklearn_crfsuite

# from sklearn_crfsuite import metrics
# from sklearn.model_selection import cross_val_predict

file = open("./assets/mod2.txt", "rb")
data = pickle.load(file)

samples = []
labels = []

exons = data[0]
introns = data[1]

i=0
while(i < len(exons)):
    for seq in exons[i]:
        samples.append([{"sequence":seq[0]}])
        labels.append([seq[-1]])
    for seq in introns[i]:
        samples.append([{"sequence":seq[0]}])
        labels.append([seq[-1]])
    i+=1

clf = sklearn_crfsuite.CRF(algorithm='lbfgs')
clf = clf.fit(samples, labels)

# pred = cross_val_predict(clf, samples, labels, cv=10)
# print("Classification_report:")
# labels = ['Intron', 'Exon', 'Neither']
# metrics.flat_f1_score(labels, pred, average='weighted', labels=labels)
# sorted_labels = sorted(labels, key=lambda name: (name[1:], name[0]))
# print(metrics.flat_classification_report(
#     labels, pred, labels=sorted_labels, digits=2))

pickle.dump(clf, open("./results/model.sav", 'wb'))