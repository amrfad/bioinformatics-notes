import replication
import motifs
import data

dsor = data.get_dosr_motif()

best_motifs = motifs.gibbs_sampler(dsor, 15, len(dsor), 10000)

print(best_motifs)