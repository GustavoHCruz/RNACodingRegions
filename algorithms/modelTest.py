import pickle
import sklearn_crfsuite

file = open("./results/model.sav", "rb")
model = pickle.load(file)

file = open("./assets/mod2.txt", "rb")
data = pickle.load(file)
exons = data[0]
introns = data[1]

acertos = [0.0,0.0,0.0]
total = [0.0,0.0,0.0]
acuracia = [0.0,0.0,0.0]
tipos = {"Exon":0,"Intron":1,"Neither":2}

i=0
while(i < 5200):
    for seq in exons[i]:
        saida = model.predict([{'sequence':seq[0]}])
        saida = str(saida[0][0])
        correto = str(seq[-1])
        if ("Exon" == correto):
            total[tipos[correto]] += 1
            if(correto == saida):
                acertos[tipos[correto]] += 1
        if ("Neither" == correto):
            total[tipos[correto]] += 1
            if(correto == saida):
                acertos[tipos[correto]] += 1            
        # print("Sendo testado:", seq[0],"\nPosições na Cadeia:",str(seq[1]), "\nResposta Correta:", seq[2])
        # print(model.predict([{'sequence':seq[0]}]))

    for seq in introns[i]:
        saida = model.predict([{'sequence':seq[0]}])
        saida = str(saida[0][0])
        correto = str(seq[-1])
        if ("Intron" == correto):
            total[tipos[correto]] += 1
            if(correto == saida):
                acertos[tipos[correto]] += 1
        if ("Neither" == correto):
            total[tipos[correto]] += 1
            if(correto == saida):
                acertos[tipos[correto]] += 1
        # print("Sendo testado:", seq[0], "\nPosições na Cadeia:", seq[1], "\nResposta Correta:", seq[2])
        # print(model.predict([{'sequence':seq[0]}]))

    i+=1

for i in range(0,len(acuracia)):
    acuracia[i] = (acertos[i]/total[i]) * 100

print("Acurácia de Exons:",acuracia[0],"Total:",int(total[0]),"Acertos:",int(acertos[0]))
print("Acurácia de Introns:",acuracia[1],"Total:",int(total[1]),"Acertos:",int(acertos[1]))
print("Acurácia de Neithers:",acuracia[2],"Total:",int(total[2]),"Acertos:",int(acertos[2]))