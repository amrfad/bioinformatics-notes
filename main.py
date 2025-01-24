import replication
import motifs
import data

motifs_1 = [
    'TCGGGGGTTTTT',
    'CCGGTGACTTAC',
    'ACGGGGATTTTC',
    'TTGGGGACTTTT',
    'AAGGGGACTTCC',
    'TTGGGGACTTCC',
    'TCGGGGATTCAT',
    'TCGGGGATTCCT',
    'TAGGGGAACTAC',
    'TCGGGTATAACC'
]

motifs_2 = [
    'GGCGTTCAGGCA',
    'AAGAATCAGTCA',
    'CAAGGAGTTCGC',
    'CACGTCAATCAC',
    'CAATAATATTCG'
]

# print(motifs.greedy_motif_search(motifs_2,3,5))

# print()

profile = {'A':  [0.4,  0.3,  0.0,  0.1,  0.0,  0.9],

'C':  [0.2,  0.3,  0.0,  0.4,  0.0,  0.1],

'G':  [0.1,  0.3,  1.0,  0.1,  0.5,  0.0],

'T':  [0.3,  0.1,  0.0,  0.4,  0.5,  0.0]}

print(motifs.pr('TCGGTA', profile))