import numpy as np

def m_count(motifs: list):
    row = len(motifs)
    col = len(motifs[0])
    count = {i: [0 for _ in range(col)] for i in 'ACGT'}
    for i in range(col):
        for j in range(row):
            try:
                symbol = motifs[j][i]
                count[symbol][i] += 1
            except KeyError:
                raise ValueError('Harusnya hanya A, C, G, T yang ada pada sekuens genom!')
    return count

def profile(motifs: list):
    count = m_count(motifs)
    col = len(count['A'])
    total = sum([count[i][0] for i in 'ACGT'])
    profile = {i:[(count[i][j] / total) for j in range(col)] for i in 'ACGT'}
    return profile

def m_consensus(motifs: list):
    consensus = ''
    count = m_count(motifs)
    col = len(count['A'])
    for t in range(col):
        max_val = -1
        sym = ''
        for i in 'ACGT':
            if count[i][t] > max_val:
                max_val = count[i][t]
                sym = i
        consensus += sym
    return consensus

def score(motifs: list):
    consensus = m_consensus(motifs)
    count = m_count(motifs)
    score = 0
    col = len(consensus)
    for i in range(col):
        for sym in 'ACGT':
            if sym != consensus[i]:
                score += count[sym][i]
    return score

def pr(text: str, profile: dict):
    pr = 1
    for i, sym in enumerate(text):
        pr *= profile[sym][i]
    return pr

def profile_most_probable_kmer(text: str, k: int, profile: list):
    n = len(text)
    prob = -1
    kmer = text[0:k]
    for i in range(n-k+1):
        pattern = text[i:i+k]
        if pr(pattern, profile) > prob:
            prob = pr(pattern, profile)
            kmer = pattern
    return kmer

def greedy_motif_search(dna: list, k: int, t: int):
    best_motifs = []
    for i in range(t):
        best_motifs.append(dna[i][0:k])
    n = len(dna[0])
    for i in range(n-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range(1, t):
            p = profile(motifs[0:j])
            motifs.append(profile_most_probable_kmer(dna[j], k, p))
        if score(best_motifs) > score(motifs):
            best_motifs = motifs
    return best_motifs

# def motifs_entropy(motifs: list):


def total_entropy(motifs: list):
    profile = profile(motifs)
    sum_h = 0
    n = len(profile['A'])
    for i in range(n):
        h = 0
        for j in 'ACGT':
            if profile[j][i] != 0:
                h += profile[j][i] * np.log2(profile[j][i])
        h = (-1) * h
        sum_h += h
    return sum_h