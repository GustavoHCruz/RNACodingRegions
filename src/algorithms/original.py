# É necessário ter instalado o biopython, o sklearn e o crfsuite
# Estes podem ser instalados pelo pip.
# Para instalar o pip use o comando sudo apt-get install python3-pip
# Para instalar o biopython use o comando pip3 install -U biopython
# Para instalar o sklearn use o comando pip3 install -U scikit-learn
# Para instalar o crfsuite use o comando pip3 install -U sklearn-crfsuite

import re
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import sklearn_crfsuite
from sklearn_crfsuite import metrics
from sklearn.model_selection import cross_val_predict
from sklearn.datasets import make_classification
from sklearn.metrics import classification_report
from sklearn import model_selection
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import multilabel_confusion_matrix
import scipy.stats
from sklearn.metrics import make_scorer
from collections import Counter
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import MultiLabelBinarizer
import numpy as num
import pickle

# recebe a lista com as posições dos éxons e faz a separação entre início, fim e exon
# Exemplo: [[0, 8, 'exon'], [106, 137, 'exon'], [198, 231, 'exon']]
def regFilter(input):
    regex = r"(\d+)\:[\>]?(\d+)"
    matches = re.finditer(regex, input, re.MULTILINE | re.IGNORECASE)

    result = []

    for matchNum, match in enumerate(matches, start=1):
        slice = []
        slice.append(int(match.group(1)))
        slice.append(int(match.group(2)))
        slice.append("exon")
        result.append(slice)
    # print(result)

    return result

# Faz a formatação para o CRF. Para cada sequencia recebida extrai as features
# sequence, minus2 e plus2. Após isso pegas as informações de sequence, minus 2 e plus2
# da sequência anterior e posterior. Se for a primeira sequência, não há nada anterior
# a ela, logo é colocado em seu lugar BOS (Begin os Sequence). Se for a últmia sequência
# não há nada posterior, então é colocado EOS (End of Sequence)
def sequence2features(sent, i):
    sequence = sent[i][0]
    minus2 = sent[i][1]
    plus2 = sent[i][2]
    features = {
        'sequence': sequence,
        'minus2': minus2,
        'plus2': plus2
    }
    if i > 0:
        sequence1 = sent[i-1][0]
        minus21 = sent[i-1][1]
        plus21 = sent[i-1][2]
        features.update({
            '-1:sequence': sequence1,
            '-1:minus2': minus21,
            '-1:plus2': plus21,
        })

    else:
        features['BOS'] = True

    if i < len(sent)-1:
        sequence1 = sent[i+1][0]
        minus21 = sent[i+1][1]
        plus21 = sent[i+1][2]
        features.update({
            '+1:sequence': sequence1,
            '+1:minus2': minus21,
            '+1:plus2': plus21
        })

    else:
        features['EOS'] = True

    return features

# Envia cada um dos possíveis introns para a função sequence2features, para serem
# formatados para o CRF
def sent2features(sent):
    return [sequence2features(sent, i) for i in range(len(sent))]


def sent2labels(sent):
    return [label for token, minus2, plus2, label in sent]


def sent2tokens(sent):
    return [token for token, minus2, plus2, label in sent]

# Nome do arquivo a ser lido com as sequencias
gb_file = "../assets/colletotrichum.gb"

# Contador de sequências aceitas do arquivo
c0 = 0

# Vetor que inclui as tags extraídas da sequência
listaTags = []

# Contador da quantidade de introns identificados
nIntron = 0

# Contador da quantidade de não introns identificados
nNotIntron = 0

# O arquivo é aberto pelo biopython e é feita uma análise para cada sequeência 
# presente no arquivo
for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
    for feature in gb_record.features:
        # Na sequência analisada, verifica se na feature seq contém alguma das seguintes letras: N, R, Y, S, W, K, M, B, D, H, V
        if(gb_record.seq.find('N')==-1 and gb_record.seq.find('R')==-1 and gb_record.seq.find('Y')==-1 and gb_record.seq.find('S')==-1 and gb_record.seq.find('W')==-1 and gb_record.seq.find('K')==-1 and gb_record.seq.find('M')==-1 and gb_record.seq.find('B')==-1 and gb_record.seq.find('D')==-1 and gb_record.seq.find('H')==-1 and gb_record.seq.find('V')==-1):
            # Verifica se a feture da sequencia possui strand diferente de -1. O -1 indica o complemento reverso
            if(feature.location.strand != -1):
                # Verifica se a sequencia possui uma feture chamada CDS. Essa feature que contém as posições já identificadas dos éxons
                if (feature.type == 'CDS'):
                    # Aumaneta em 1 a quantidade total de sequências aceitas
                    c0 += 1
                    # Chama a função regFilter para obter as posições dos éxons já formatadas em um vetor
                    positions = regFilter(str(feature.location))
                    # Se o valor na posição [0][0] for > 0, significa que existe um intron da posição
                    # 0 até a posição position[0][0], logo é incluído no vetor de posições essa informção do intron
                    if(positions[0][0] > 0):
                        positions.insert(0, [0, positions[0][0], "intron"])
                    # Se o valor na posição final [len(positions)-1][1] for < que o tamanho total da sequência, 
                    # significa que existe um intron na posição [len(positions)-1][1] até o tamanho final da sequência,
                    # logo é incluído no vetor de posições essa informção do intron
                    if(positions[len(positions)-1][1] < len(gb_record.seq)):
                        positions.append(
                            [positions[len(positions)-1][1], len(gb_record.seq), "intron"])
                    # Obtém o tamanho total do vetor se sequências após a inserção dos íntrons
                    total = len(positions)
                    cont = 0
                    # É percorrido o vetor de posições, para que os introns sejam
                    # inseridos entre dois éxons seguidos
                    while(cont != total-1):
                        c1 = positions[cont+1][0]
                        c2 = positions[cont][1]
                        if(positions[cont+1][0] != positions[cont][1]):
                            positions.insert(cont+1, [c2, c1, "intron"])
                            total += 1
                        cont += 1
                    
                    # Obtém a sequência completa
                    content = str(gb_record.seq)
                    # Se o vetor positions ainda possuir algum valor, é removido o
                    # valor que se encontra na última posição. Dessa forma evita-se
                    # analisar trechos que possam estar cortados
                    if positions:
                        positions.pop(len(positions)-1)
                    # Se o vetor positions ainda possuir algum valor, é removido o
                    # valor que se encontra na primeira posição. Dessa forma evita-se
                    # analisar trechos que possam estar cortados
                    if positions:
                        positions.pop(0)
                    
                    # Inicializa listaAux com um vetor vazio
                    listaAux = []
                    # Verifica se positions ainda possui algum valor
                    if positions:
                        # Obtém o valor da última posição da sequência
                        tamSeq = len(gb_record.seq) - 1
                        listaSeqs = []
                        posIntr = []
                        # Cria um vetor chamado content2 contendo as posições que iniciam com GT
                        content2 = [m.start() for m in re.finditer('GT', content)]
                        # Cria um vetor chamado content2 contendo as posições que iniciam com GT
                        content3 = [m.start() for m in re.finditer('AG', content)]
                        possibles = []
                        rt = False
                        # Analisa cada valor do vetor content2
                        for n in content2:
                            # Analisa cada valor do vetor content3
                            for o in content3:
                                rt = False
                                # Se a posição de GT vier antes de AG
                                if n < o+2:
                                    # Analisa cada valor do vetor positions
                                    for pos in positions:
                                        # Verifica se a posição corresponde a um intron
                                        if pos[2] == 'intron':
                                            # Verifica se a posição inicial do intron é igual a posição de GT analisada
                                            # E se a posição final do intron é igual a posição de AG analisada + 2
                                            if pos[0] == n and pos[1] == o+2:
                                                # Verifica se a posição de GT - 2 é maior ou igual a zero
                                                if(n-2 >= 0):
                                                    # Se for é adicionado em p1 os dois nucleotídeos anteriores a GT da posição n
                                                    p1 = content[n-2:n]
                                                else:
                                                    # Se for menor que zero é adicionada a string BEGIN
                                                    p1 = 'BEGIN'
                                                # Verifica se a posição de AG + 4 é menor ou igual a posição final da sequência
                                                if(o+4 <= len(gb_record.seq)-1):
                                                    # Se for é adicionado em p2 os dois nucleotídeos posteriores a AG da posição o
                                                    p2 = content[o+2:o+4]
                                                else:
                                                    # Se for maior que o final da sequência é adicionada a string END
                                                    p2 = 'END'
                                                # Adiciona em listaAux o conteúdo da primeira sequência analisada, que contém
                                                # a sequência completa, os valores de p1, p2 e INTRON
                                                listaAux.append(
                                                    (content[n:o+2], p1, p2, 'INTRON'))
                                                rt = True
                                                nIntron += 1
                                    # Se rt for false é repetido o processo anterior mas substituindo o INTROn por NOT-INTRON
                                    if(rt == False):
                                        if(n-2 >= 0):
                                            p1 = content[n-2:n]
                                        else:
                                            p1 = 'BEGIN'
                                        if(o+4 <= len(gb_record.seq)-1):
                                            p2 = content[o+2:o+4]
                                        else:
                                            p2 = 'END'
                                        listaAux.append(
                                            (content[n:o+2], p1, p2, 'NOT-INTRON'))
                                        nNotIntron += 1
                    listaTags.append(listaAux)

listaAux = []
print('Size 1: ', c0)

X_train = [sent2features(s) for s in listaTags]
y_train = [sent2labels(s) for s in listaTags]

print("CRF 2.0")
clf = sklearn_crfsuite.CRF(algorithm='lbfgs', c1=0.07, c2=0.09,
                           max_iterations=100, all_possible_transitions=True)
# clf = clf.fit(X_train, y_train)
pred = cross_val_predict(clf, X_train, y_train, cv=10)
# pred = clf.predict(X_train)
print("Classification_report:")
labels = ['NOT-INTRON', 'INTRON']
metrics.flat_f1_score(y_train, pred, average='weighted', labels=labels)
sorted_labels = sorted(labels, key=lambda name: (name[1:], name[0]))
print(metrics.flat_classification_report(
    y_train, pred, labels=sorted_labels, digits=2))

# O trecho abaixo utiliza a biblioteca pickle para salvar o resultado do treinamento em um arquivo chamado finalized_model.sav
# num.savetxt('pred.txt', pred, delimiter=" ", fmt="%s")
# filename = 'finalized_model.sav'
# pickle.dump(clf, open(filename, 'wb'))

print('Quantidade de Introns:', nIntron)
print('Quantidade de não introns:', nNotIntron)