# Write a program to identify restriction recognition sites in a given DNA sequence. [Hint:
# Take a small sequence from the list of sequences available in ReBase Database.]

complement = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

NUCLEOTIDES = [ 'A', 'T', 'C', 'G' ]

forward_dna_strand = input("Enter Forward DNA strand: ")

for nucleotides in forward_dna_strand:
    if nucleotides.upper() not in NUCLEOTIDES:
        print("Unknown DNA sequence")
        quit()

forward_dna_strand = [ i.upper() for i in forward_dna_strand ]
complement_dna_strand = [ complement[i] for i in forward_dna_strand ]
reversed_complement_strand = complement_dna_strand[::-1]
print("forward_dna_strand : ","".join(forward_dna_strand))
print("complement_dna_strand : ","".join(complement_dna_strand))
print("reversed_complement_strand : ","".join(reversed_complement_strand))


dna_length = len(forward_dna_strand)
restriction_recognition_site = []
restriction_recognition_site_index = []

restriction_recognition_site_size = [4,5,6,7,8]
for sz in restriction_recognition_site_size:
    for i in range(dna_length-sz+1):
        if forward_dna_strand[i:i+sz]!=reversed_complement_strand[dna_length-sz-i:dna_length-i]:
            continue
        else:
            restriction_recognition_site.append(forward_dna_strand[i:i+sz-1])
            restriction_recognition_site_index.append((i,i+sz))

final_output = []
for site in restriction_recognition_site:
    final_output.append("".join(site))

print("Total Pallindromic sites : ",len(final_output))
print("\nS.No.     Position        Site")
for i in range(len(final_output)):
    print(f"{i+1}         {restriction_recognition_site_index[i]}          {final_output[i]}")