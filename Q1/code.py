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

# i.e.
# gtttcattataccagtttagatctatcgacagggcgttgagtgtgtgcttactcacggctggcatgtaggtaacagtagtggggaagcgtaacatctgaggcctgactcacatatagagtgtcgaccaaggggtgaagcatcatacgccatacaggcccctagcgaaacgacctagtctaaagacacacgagaatgaaacccgtggacttggttacagcgtaataatctggtcagagctggtccggcgctggcgatgtaccttacgccactgcaaaccggctttgcagagaacatctgggtacattcccgtgtcatgtcaaagcaggtgattcccgcgaaaaacaattaacgacgcatttgctattgacgaagtcctagttctccgaattgagcgggagacatatgatgtcgagactgcaggaaccgaattatcctgtccgcagatccaatagctcacagaggtaaggggagtgtgatggtgccctagggtgtttgaacg

AMINO_ACID = {
    "UUU":"Phe","UUC":"Phe","UUA":"Leu","UUG":"Leu",
    "UGU":"Cys","UGC":"Cys","UGA":"-","UGG":"Trp",
    "UAU":"Tyr","UAC":"Tyr","UAA":"-","UAG":"-",
    "UCU":"Ser","UCC":"Ser","UCA":"Ser","UCG":"Ser",
    "CUU":"Leu","CUC":"Leu","CUG":"Leu","CUA":"Leu",
    "CGU":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "CCU":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "CAU":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "GUU":"Val","GUC":"Val","GUA":"Val","GUG":"Val",
    "GCU":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "GGU":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
    "GAU":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "AUU":"Ile","AUC":"Ile","AUA":"Ile","AUG":"Met",
    "ACU":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "AGU":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "AAU":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
}

RNA = { # assuming forward DNA strand
    'A':'A',
    'T':'U',
    'C':'C',
    'G':'G'
}

NUCLEOTIDES = [ 'A', 'T', 'C', 'G' ]

def verify_dna_strand(dna_strand):
    '''
    This function will verify dna_strand.
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


def translate(rna_strand):
    '''
    This function convert given rna strand string into its protein strand string
    '''
    MET_INDEX = rna_strand.find('AUG')
    if MET_INDEX == -1:
        return "Start Codon Not Found"

    protein_strand = ""
    t_rna_strand = rna_strand[MET_INDEX:]
    t_rna_length = len(t_rna_strand)
    amino_acid_size = 3
    mod = t_rna_length%amino_acid_size

    if mod==0:
        for i in range(0,t_rna_length,amino_acid_size):
            codon = t_rna_strand[i:i+amino_acid_size].upper()
            protein = AMINO_ACID[codon]
            protein_strand += protein
            if protein=='-':
                return protein_strand
    
    elif mod==1:
        for i in range(0,t_rna_length-1,amino_acid_size):
            codon = t_rna_strand[i:i+amino_acid_size].upper()
            protein = AMINO_ACID[codon]
            protein_strand += protein
            if protein=='-':
                return protein_strand
        protein_strand+=t_rna_strand[-1]
    
    else:
        for i in range(0,t_rna_length-2,amino_acid_size):
            codon = t_rna_strand[i:i+amino_acid_size].upper()
            protein = AMINO_ACID[codon]
            protein_strand += protein
            if protein=='-':
                return protein_strand
        protein_strand+=t_rna_strand[-2:]
    
    return protein_strand


if __name__=="__main__":
    dna_strand = input()
    if not verify_dna_strand(dna_strand):
        print("Unknown DNA Sequence")
        quit()
    rna_strand = transcript(dna_strand)
    protein_strand = translate(rna_strand)
    print("RNA strand:")
    print(rna_strand)
    print("PROTEIN strand:")
    print(protein_strand)