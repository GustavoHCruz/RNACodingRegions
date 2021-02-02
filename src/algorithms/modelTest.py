# Main ===================================================================================================================
# Libraries and importations ===========================
import pickle               # Used for file operations

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


# Cross-validate prediction. Results are printed on the console
pred = clf.predict(testSamples)
predCV = cross_val_predict(clf, trainSamples, trainLabels, cv=10)
x = ['Intron', 'Exon', 'Neither']
# metrics.flat_f1_score(testLabels, pred, average='weighted', labels=x)
sorted_labels = sorted(x, key=lambda name: (name[1:], name[0]))
print("Classification Report:")
print(metrics.flat_classification_report(testLabels, pred, labels=sorted_labels, digits=2))
print("Classification Report CV:")
print(metrics.flat_classification_report(trainLabels, predCV, labels=sorted_labels, digits=2))
# ========================================================================================================================