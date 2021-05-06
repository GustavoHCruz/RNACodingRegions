from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

result = NCBIWWW.qblast("blastp", "nt", "SIPRPQSDPFRRSRSAVVGLRKGHVDSDNTSTSIVGRPRHHGIMIGMG")

blast_record = NCBIXML.read(result)

E_VALUE_THRESH = 0.04

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")