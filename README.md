# Descrição
Repositório contendo os algoritmos desenvolvidos no Projeto de Iniciação Científica. O PIC durou um ano e teve como objetivo criar um modelo de aprendizado de máquina (AM) capaz de prever regiões de _splicing_ no RNA. Foram usados registros dos fungos _Colletotrichum_ e _Diaporthe_ para a criação do modelo. Detalhes sobre o desenvolvimento do projeto, resultados e informações adicionais estão disponíveis no arquivo **Relatório.pdf**.

# Bibliotecas necessárias
Para realizar um novo treinamento, é necessário ter as seguintes bibliotecas instaladas:
- [Biopython] (https://github.com/biopython/biopython): pip install biopython
- [Sklearn] (https://github.com/scikit-learn/scikit-learn): pip install scikit-learn
- [Sklearn-crfsuite] (https://github.com/TeamHG-Memex/sklearn-crfsuite/blob/master/docs/index.rst): pip install sklearn-crfsuite
- 
# Como usar
O treinamento é feito por partes.
- O primeiro passo é executar o módulo 1 (**Module 1 - Selection.py**), onde se deve dizer qual o nome do arquivo de entrada (arquivo provindo do [GenBank] (https://www.ncbi.nlm.nih.gov/genbank/)) e o nome do arquivo de saída (ambas as informações devem ser alteradas nas linhas **55** e **56** do código).
- O segundo passo é executar o módulo 2 (**Module 2 - Preparation.py**), informando o nome do arquivo de dados gerado no módulo anterior e um nome para o arquivo de saída deste módulo (ambas as informações devem ser alteradas nas linhas **23** e **24** do código).
- O terceiro passo é executar o módulo 3 (**Module 3 - Training.py**), informando o nome do arquivo de dados gerado no módulo 2, e um nome para o arquivo de saída deste módulo e para o arquivo onde o modelo treinado será armazenado (todas as informações devem ser alteradas nas linhas **9**, **10** e **11** do código).
- Com isso, pode-se executar o módulo 4 (**Module 4 - Application.py**), informando o nome do arquivo do modelo treinado e a sequência a ser verificada pelo modelo. O resultado é impresso no terminal.
Existem ainda módulos para a avaliar a precisão geral do modelo, como o **Test - CrossValidation.py** (precisa do nome do arquivo do modelo), que usa o método de _crossvalidation_ para gerar um resultado da precisão do modelo e o **Test - Final Sequence Comparation.py** (precisa do nome do arquivo do modelo e do nome do arquivo do módulo 1), que avaliará quantas sequências o modelo prediz corretamente tendo como base os dados gerados pelo módulo 1.

# Autores
- Gustavo Henrique Ferreira Cruz
- Vinícius Menossi
