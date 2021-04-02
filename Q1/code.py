# Write a code to accept any DNA sequence (of varying lengths) and produce as output the
# corresponding RNA strand synthesized and protein strand synthesized?

# Sample DNA sequence:
# gtttcattataccagtttagatctatcgacagggcgttgagtgtgtgcttactcacggct
# ggcatgtaggtaacagtagtggggaagcgtaacatctgaggcctgactcacatatagagt
# gtcgaccaaggggtgaagcatcatacgccatacaggcccctagcgaaacgacctagtcta
# aagacacacgagaatgaaacccgtggacttggttacagcgtaataatctggtcagagctg
# gtccggcgctggcgatgtaccttacgccactgcaaaccggctttgcagagaacatctggg
# tacattcccgtgtcatgtcaaagcaggtgattcccgcgaaaaacaattaacgacgcattt
# gctattgacgaagtcctagttctccgaattgagcgggagacatatgatgtcgagactgca
# ggaaccgaattatcctgtccgcagatccaatagctcacagaggtaaggggagtgtgatgg
# tgccctagggtgtttgaacg


# We can directly obtain protein_strand from dna_strand, no need for converting it to
# rna_strand as only changes in the triplet will be that 'T' is replaced by 'U'


AMINO_ACID = {
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}


RNA = {
    'A':'U',
    'T':'A',
    'C':'G',
    'G':'C'
}

NUCLEOTIDES = [ 'A', 'T', 'C', 'G' ]

def verify_dna_strand(dna_strand):
    '''
    '''
    for i in range(len(dna_strand)):
        if dna_strand[i].upper() not in NUCLEOTIDES:
            return False
    return True


def transcript(dna_strand):
    '''
    This function convert given dna strand string into rna strand string
    '''
    rna_strand = ""
    for i in range(len(dna_strand)):
        rna_strand += RNA[dna_strand[i].upper()]
    return rna_strand


def translate(dna_strand):
    '''
    This function convert given dna strand string into its protein strand string
    '''

if __name__=="__main__":
    dna_strand = input()
    if not verify_dna_strand(dna_strand):
        print("Unknown DNA Sequence")
        quit()
    rna_strand = transcript(dna_strand)
    protein_strand = translate(dna_strand)
    print("RNA strand:")
    print(rna_strand)
    print("PROTEIN strand:")
    print(protein_strand)