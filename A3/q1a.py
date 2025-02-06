from Bio import SeqIO
import pandas as pd

# load genome sequence
genome_file = "C:\\Users\\tinas\\Downloads\\Vibrio_cholerae.GFC_11.dna.toplevel.fa"
genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# load GFF3 gene annotation
gff_file = "C:\\Users\\tinas\\Downloads\\Vibrio_cholerae.GFC_11.37.gff3"
annotations = pd.read_csv(gff_file, sep='\t', comment='#', header=None)

# filter the GFF3 annotations to get CDS features on the forward strand
cds_annotations = annotations[annotations[2] == 'CDS']
cds_annotations = cds_annotations[cds_annotations[6] == '+']  # forward strand only

# i) calculate average length of intergenic regions
intergenic_lengths = []
total_non_coding_length = 0  # stores the total length of non-coding regions

# check for region before the first CDS
first_cds_start = cds_annotations.iloc[0, 3]
if first_cds_start > 0:
    intergenic_lengths.append(first_cds_start)  # region from 0 to the start of the first CDS
    total_non_coding_length += first_cds_start  # add the length of this initial intergenic region

# process the intergenic regions between consecutive CDS features
for i in range(1, len(cds_annotations)):
    previous_end = cds_annotations.iloc[i - 1, 4]
    current_start = cds_annotations.iloc[i, 3]
    
    # intergenic region: between the end of the previous CDS and start of the next CDS
    if current_start > previous_end:
        intergenic_lengths.append(current_start - previous_end)
        total_non_coding_length += (current_start - previous_end)  # add the length of this intergenic region to the total

print(f"Total length of intergenic regions: {total_non_coding_length} bp")
print(f"Number of intergenic regions: {len(intergenic_lengths)}")
avg_intergenic_length = sum(intergenic_lengths) / len(intergenic_lengths) if intergenic_lengths else 0
print(f"Average length of intergenic regions: {avg_intergenic_length} bp")

# ii) calculate average length of CDS
cds_lengths = cds_annotations[4] - cds_annotations[3]
print(f"Total length of CDS: {cds_lengths.sum()}")
print(f"Number of CDS annotations: {len(cds_annotations)}")
print(f"Average length of CDS: {cds_lengths.mean()} bp")

# iii) calculate the nucleotide frequency in intergenic regions
# collect intergenic sequences
intergenic_sequences = []
for i in range(1, len(cds_annotations)):
    previous_end = cds_annotations.iloc[i - 1, 4]
    current_start = cds_annotations.iloc[i, 3]
    
    if current_start > previous_end:
        contig = cds_annotations.iloc[i, 0]  # Contig/Chromosome
        intergenic_seq = genome[contig][previous_end:current_start].seq
        intergenic_sequences.append(str(intergenic_seq))

# add the sequence before the first CDS (if it exists)
if first_cds_start > 0:
    contig = cds_annotations.iloc[0, 0]  # Contig/Chromosome
    intergenic_seq = genome[contig][0:first_cds_start].seq
    intergenic_sequences.insert(0, str(intergenic_seq))

# join all intergenic sequences into one long string
intergenic_sequence = ''.join(intergenic_sequences)

# count nucleotide frequencies
nucleotide_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
total_nucleotides = len(intergenic_sequence)

# count each nucleotide in the intergenic sequence
for nt in intergenic_sequence:
    if nt in nucleotide_counts:
        nucleotide_counts[nt] += 1

nucleotide_frequencies = {nt: count / total_nucleotides for nt, count in nucleotide_counts.items()}
print("Nucleotide frequency table for intergenic regions:")
for nt, freq in nucleotide_frequencies.items():
    print(f"{nt}: {freq}")

# iv) calculate the codon frequency in CDS
start_codons = ['ATG', 'GTG', 'TTG']
stop_codons = ['TAA', 'TAG', 'TGA']

# create dictionaries for counting start, stop, and internal codons
start_codon_counts = {codon: 0 for codon in start_codons}
stop_codon_counts = {codon: 0 for codon in stop_codons}
internal_codon_counts = {}

# go through each CDS
for index, row in cds_annotations.iterrows():
    contig = row[0]  # Chromosome/contig
    start = row[3] - 1
    end = row[4]

    # grab the CDS sequence
    cds_seq = genome[contig][start:end].seq.upper()

    # separate the sequence into codons
    codons = [str(cds_seq[i:i+3]) for i in range(0, len(cds_seq), 3)]
    
    # skip incomplete last codon if it exists
    if len(cds_seq) % 3 != 0:
        codons = codons[:-1]
    
    # find and count start codon
    if codons[0] in start_codons:
        start_codon_counts[codons[0]] += 1
    
    # find and count stop codon
    if codons[-1] in stop_codons:
        stop_codon_counts[codons[-1]] += 1
    
    # count internal codons (excluding start and stop)
    for codon in codons[1:-1]:  # exclude the first and last codon
        if codon not in internal_codon_counts:
            internal_codon_counts[codon] = 0
        internal_codon_counts[codon] += 1

# compute total counts for frequencies
total_start_codons = sum(start_codon_counts.values())
total_stop_codons = sum(stop_codon_counts.values())
total_internal_codons = sum(internal_codon_counts.values())

# compute frequency tables
start_codon_frequencies = {codon: count / total_start_codons for codon, count in start_codon_counts.items()}
stop_codon_frequencies = {codon: count / total_stop_codons for codon, count in stop_codon_counts.items()}
internal_codon_frequencies = {codon: count / total_internal_codons for codon, count in internal_codon_counts.items()}

print("START Codon Frequency Table:")
for codon, freq in start_codon_frequencies.items():
    print(f"{codon}: {freq}")

print("\nSTOP Codon Frequency Table:")
for codon, freq in stop_codon_frequencies.items():
    print(f"{codon}: {freq}")

print("\nInternal Codon Frequency Table:")
for codon, freq in internal_codon_frequencies.items():
    print(f"{codon}: {freq}")
