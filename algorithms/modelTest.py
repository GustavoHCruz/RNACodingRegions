import pickle
import sklearn_crfsuite

file = open("./results/model.sav", "rb")
model = pickle.load(file)

file = open("./assets/mod2.txt", "rb")
data = pickle.load(file)
exons = data[0]
introns = data[1]

hits = [0.0,0.0,0.0]
total = [0.0,0.0,0.0]
accuracy = [0.0,0.0,0.0]
types = {"Exon":0,"Intron":1,"Neither":2}

i=0
while(i < 5200):
    for seq in exons[i]:
        output = model.predict([{'sequence':seq[0]}])
        output = str(output[0][0])
        correctOne = str(seq[-1])
        if ("Exon" == correctOne):
            total[types[correctOne]] += 1
            if(correctOne == output):
                hits[types[correctOne]] += 1
        if ("Neither" == correctOne):
            total[types[correctOne]] += 1
            if(correctOne == output):
                hits[types[correctOne]] += 1            
        # print("Being Tested:", seq[0],"\nSequence Positions:",str(seq[1]), "\nCorrectly Answer:", seq[2])
        # print(model.predict([{'sequence':seq[0]}]))

    for seq in introns[i]:
        output = model.predict([{'sequence':seq[0]}])
        output = str(output[0][0])
        correctOne = str(seq[-1])
        if ("Intron" == correctOne):
            total[types[correctOne]] += 1
            if(correctOne == output):
                hits[types[correctOne]] += 1
        if ("Neither" == correctOne):
            total[types[correctOne]] += 1
            if(correctOne == output):
                hits[types[correctOne]] += 1
        # print("Being Tested:", seq[0], "\nSequence Positions:", seq[1], "\nCorrectly Answer:", seq[2])
        # print(model.predict([{'sequence':seq[0]}]))

    i+=1

for i in range(0,len(accuracy)):
    accuracy[i] = (hits[i]/total[i]) * 100

print("Exons accuracy:",accuracy[0],"Total:",int(total[0]),"hits:",int(hits[0]))
print("Introns accuracy:",accuracy[1],"Total:",int(total[1]),"hits:",int(hits[1]))
print("Neithers accuracy:",accuracy[2],"Total:",int(total[2]),"hits:",int(hits[2]))