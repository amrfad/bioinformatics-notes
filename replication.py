import matplotlib.pyplot as plt

def pattern_count(text: str, pattern: str):
    """
    Menghitung berapa kali suatu pola (pattern) muncul dalam teks (text).

    Args:
        text (str): Teks atau urutan yang akan dicari polanya.
        pattern (str): Pola yang ingin dihitung kemunculannya.

    Returns:
        int: Jumlah kemunculan pola dalam teks.
    """
    count = 0
    len_text = len(text)
    len_pattern = len(pattern)
    for i in range(len_text-len_pattern+1):
        if text[i:i+len_pattern] == pattern:
            count = count+1
    return count

def frequency_map(text: str, k: str):
    """
    Membuat peta frekuensi (frequency map) untuk semua kemungkinan k-mer (substring panjang k) dalam teks.

    Args:
        text (str): Teks atau urutan yang akan dianalisis.
        k (int): Panjang k-mer yang ingin dicari.

    Returns:
        dict: Dictionary dengan kunci berupa k-mer dan nilai berupa frekuensi kemunculannya.
    """
    freq = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:i+k]
        freq[pattern] = 0
    for kmer in freq:
        freq[kmer] = pattern_count(text, kmer)
    return freq

def frequent_words(text: str, k: int):
    """
    Mencari semua k-mer yang paling sering muncul dalam teks.

    Args:
        text (str): Teks atau urutan yang akan dianalisis.
        k (int): Panjang k-mer yang ingin dicari.

    Returns:
        list: Daftar k-mer yang paling sering muncul.
    """
    freq = frequency_map(text, k)
    words = []
    m = max(freq.values())
    for kmer in freq:
        if m == freq[kmer]:
            words.append(kmer)
    return words

def reverse(pattern: str):
    """
    Membalikkan urutan suatu pola (pattern).

    Args:
        pattern (str): Pola yang ingin dibalik.

    Returns:
        str: Pola yang sudah dibalik.
    """
    return pattern[::-1]

def complement(pattern: str):
    """
    Mencari komplemen dari suatu pola DNA (A ↔ T, C ↔ G).

    Args:
        pattern (str): Pola DNA yang ingin dicari komplemennya.

    Returns:
        str: Pola DNA yang sudah dikomplemen.
    """
    complement = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    pattern = list(pattern)
    for i, base in enumerate(pattern):
        pattern[i] = complement[base]
    pattern = ''.join(pattern)
    return pattern

def reverse_complement(pattern: str):
    """
    Mencari reverse complement dari suatu pola DNA (membalikkan dan mengkomplemen).

    Args:
        pattern (str): Pola DNA yang ingin dicari reverse complement-nya.

    Returns:
        str: Pola DNA yang sudah dibalik dan dikomplemen.
    """
    return complement(reverse(pattern))

def pattern_matching(genome: str, pattern: str):
    """
    Mencari semua posisi awal kemunculan suatu pola dalam genom.

    Args:
        genome (str): Genom atau urutan DNA yang akan dicari.
        pattern (str): Pola yang ingin dicari posisinya.

    Returns:
        list: Daftar posisi awal kemunculan pola dalam genom.
    """
    position = []
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            position.append(i)
    return position

def symbol_array(genome: str, symbol: str):
    """
    Menghitung jumlah kemunculan suatu simbol dalam setengah genom (half-genome) 
    untuk setiap posisi dalam genom.

    Args:
        genome (str): Genom atau urutan DNA yang akan dianalisis.
        symbol (str): Simbol yang ingin dihitung kemunculannya.

    Returns:
        list: Daftar jumlah kemunculan simbol dalam setengah genom untuk setiap posisi.
    """
    n = len(genome)
    extended_genome = genome + genome[0:n//2]
    array = []
    array.append(pattern_count(extended_genome[0:(n//2)], symbol))
    for i in range(1, n):
        count = array[i-1]
        if extended_genome[i-1] == symbol: count -= 1
        if extended_genome[i+(n//2)-1] == symbol: count += 1
        array.append(count)
    return array

def skew_array(genome: str):
    """
    Menghitung skew (perbedaan antara jumlah G dan C) untuk setiap posisi dalam genom.

    Args:
        genome (str): Genom atau urutan DNA yang akan dianalisis.

    Returns:
        list: Daftar nilai skew untuk setiap posisi dalam genom.
    """
    skew = [0]
    n = len(genome)
    for i in range(n):
        if genome[i] == 'C':
            skew.append(skew[i] - 1)
        elif genome[i] == 'G':
            skew.append(skew[i] + 1)
        else:
            skew.append(skew[i])
    return skew

def minimum_skew(genome: str):
    """
    Mencari posisi dengan nilai skew minimum dalam genom.

    Args:
        genome (str): Genom atau urutan DNA yang akan dianalisis.

    Returns:
        list: Daftar posisi dengan nilai skew minimum.
    """
    skew = skew_array(genome)
    mins = min(skew)
    n = len(skew)
    position = []
    for i in range(n):
        if skew[i] == mins: position.append(i)
    return position

def hamming_distance(str1: str, str2: str):
    """
    Menghitung jarak Hamming antara dua string (jumlah posisi di mana karakter berbeda).

    Args:
        str1 (str): String pertama.
        str2 (str): String kedua.

    Returns:
        int: Jarak Hamming antara dua string.

    Raises:
        ValueError: Jika panjang kedua string tidak sama.
    """
    if len(str1) != len(str2):
        raise ValueError('Panjang kedua string harus sama!')
    n = len(str1)
    distance = 0
    for i in range(n):
        if str1[i] != str2[i]: distance += 1
    return distance

def approx_pattern_matching(genome: str, pattern: str, d: int):
    """
    Mencari semua posisi awal kemunculan suatu pola dalam genom dengan toleransi 
    perbedaan (mismatch) sebanyak d.

    Args:
        genome (str): Genom atau urutan DNA yang akan dicari.
        pattern (str): Pola yang ingin dicari posisinya.
        d (int): Jumlah maksimum mismatch yang diizinkan.

    Returns:
        list: Daftar posisi awal kemunculan pola dalam genom dengan toleransi mismatch.
    """
    position = []
    n_genome = len(genome)
    n_pattern = len(pattern)
    for i in range(n_genome-n_pattern+1):
        if hamming_distance(genome[i:i+n_pattern], pattern) <= d:
            position.append(i)
    return position

def approx_pattern_count(genome: str, pattern: str, d: int):
    """
    Menghitung berapa kali suatu pola muncul dalam teks dengan toleransi perbedaan (mismatch) sebanyak d.

    Args:
        genome (str): Genom atau urutan DNA yang akan dicari polanya.
        pattern (str): Pola yang ingin dihitung kemunculannya.
        d (int): Jumlah maksimum mismatch yang diizinkan.

    Returns:
        int: Jumlah kemunculan pola dalam teks dengan toleransi mismatch.
    """
    count = 0
    len_text = len(genome)
    len_pattern = len(pattern)
    for i in range(len_text-len_pattern+1):
        if hamming_distance(genome[i:i+len_pattern], pattern) <= d:
            count = count+1
    return count

def visualize_symbol_array(array: list, symbol:str=None):
    """
    Membuat visualisasi scatter plot dari symbol array.

    Args:
        array (list): Array yang berisi jumlah simbol untuk setiap posisi.
        symbol (str, optional): Simbol yang divisualisasikan. Default: None.
    """
    X = [i for i in range(len(array))]
    y = array
    plt.scatter(X, y, color='blue', s=0.001)
    plt.xlabel('Posisi genome')
    if symbol:
        plt.ylabel(f'Jumlah {symbol} dalam half-genome berdasarkan posisi x')
        plt.title(f'Grafik persebaran {symbol} dalam genom')
    else:
        plt.ylabel('Jumlah Simbol')
        plt.title('Grafik Persebaran Simbol')
    plt.show()

def visualize_skew(genome: str):
    """
    Membuat visualisasi scatter plot dari skew array.

    Args:
        genome (str): Genom atau urutan DNA yang akan divisualisasikan skew-nya.
    """
    skew = skew_array(genome)
    X = [i for i in range(len(skew))]
    y = skew
    plt.scatter(X, y, color='blue', s=0.001)
    plt.xlabel('Posisi Genom')
    plt.ylabel('Skew')
    plt.title('Visualisasi Skew Diagram')
    plt.show()