# Functions ==============================================================================================================
# String -> List of lists
# Transform the entry string to a list with just the important data
# Entry example: 'join{[<0:29](+), [130:161](+), [230:>243](+)}'
# Output example: [[0,28],[130,160],[230,243]]
def make_exons_list(seq):
    seq = str(seq)
    seq = seq.split("{")
    seq = seq[-1]
    for char in "}][><)(+ ":
        seq = seq.replace(char,"")
    seq = seq.split(",")
    i = 0
    while(i < len(seq)):
        seq[i] = seq[i].split(":")
        seq[i][0] = int(seq[i][0])
        seq[i][1] = int(seq[i][1])
        seq[i][1] = seq[i][1] - 1
        i += 1
    return seq
assert(make_exons_list("join{[<12:50](+), [102:130](+), [186:>201](+)}") == [[12,49],[102,129],[186,200]])
assert(make_exons_list("join{[<0:40](+), [230:>243]}") == [[0,39],[230,242]])
assert(make_exons_list("join{[47:90](+)}") == [[47,89]])

# List of lists and integer -> list of lists
# Using the exons positions create the list of the introns positions
# Entry example: [[0,28],[130,160],[230,243]],244
# Output example: [[29,129],[161,229]]
def make_introns_list(exons_list,length):
    i = 0
    seq = []
    if(exons_list[i][0] != 0):
        seq.append([0,exons_list[i][0]-1])
    while(i < len(exons_list) - 1):
        seq.append([exons_list[i][1]+1,exons_list[i+1][0]-1])
        i += 1
    if(exons_list[-1][1] != length - 1):
        seq.append([exons_list[-1][1]+1,length-1])

    return seq
assert(make_introns_list([[0,28],[130,160],[230,243]],244) == [[29,129],[161,229]])
assert(make_introns_list([[0,99]],100) == [])
assert(make_introns_list([[20,99]],100) == [[0,19]])
assert(make_introns_list([[0,60]],100) == [[61,99]])
assert(make_introns_list([[20,60]],100) == [[0,19],[61,99]])
# ========================================================================================================================

# Main ===================================================================================================================
# Libraries and importations ===========================
from Bio import SeqIO
import pickle # Used for file operations
# ======================================================

# Genbank archive to be used (you can edit the file name)
genbank_archive = open("./assets/diaporthe.gb","r")

# Variable that contains all of the processed data
data = []

# Check every register in archive
for register in SeqIO.parse(genbank_archive, "genbank"):
    # Check if the register is valid (CDS means that if the sequence is reversed)
    if(register.features[-1].type == "CDS"):
        # Transform the sequence of the register to a String type
        seq = str(register.seq)

        # Create the exons list from the genbank register
        exons_list = make_exons_list(register.features[-1].location)

        # Create the introns list from the genbank register
        introns_list = make_introns_list(exons_list,len(seq))

        # Lists that cointains the sequences of introns and exons
        introns,exons = [],[]

        # Put the sequence in variables using his respectives indexs
        for x in exons_list:
            exons.append([seq[(x[0]):(x[1])+1]])
        for x in introns_list:
            introns.append([seq[(x[0]):(x[1])+1]])

        # Append all of the processed data in data variable
        data.append([seq,exons_list,exons,introns_list,introns])

        # Saves the result in the file (you can edit the name of destination file)
        file = open("./assets/diaporthe_mod1.txt","wb")
        pickle.dump(data,file)
# ========================================================================================================================