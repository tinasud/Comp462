from Bio import SeqIO

# Helper function to determine whether it is a match or mistmatch
def get_score(match, mismatch, Si, Tj):
    if Si == Tj:
        return match
    else:
        return mismatch

# Function which implements the slippage aware alignment algorithm
def get_slippage_aware_alignment(S, T, matchScore, mismatchScore, cs, cn):
    m = len(S)
    n = len(T)
    X = []

    # Initialize all entries in X to 0 (X is the scoring matrix)
    for _ in range(m + 1):
        X.append([0] * (n + 1))

    # Fill in first row of X
    for i in range(1, m + 1):
        if i > 1 and S[i - 1] == S[i - 2]:
            X[i][0] = X[i - 1][0] + cs
        else:
            X[i][0] = X[i - 1][0] + cn

    # Fill in first column of X
    for j in range(1, n + 1):
        if j > 1 and T[j - 1] == T[j - 2]:
            X[0][j] = X[0][j - 1] + cs
        else:
            X[0][j] = X[0][j - 1] + cn

    # Fill the rest of the table using the slippage-aware alignment scoring scheme
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            case1 = X[i - 1][j - 1] + get_score(matchScore, mismatchScore, S[i - 1], T[j - 1])
            
            if i > 1 and S[i - 1] == S[i - 2]:
                case2 = X[i - 1][j] + cs
            else:
                case2 = X[i - 1][j] + cn
            
            if j > 1 and T[j - 1] == T[j - 2]:
                case3 = X[i][j - 1] + cs
            else:
                case3 = X[i][j - 1] + cn

            X[i][j] = max(case1, case2, case3)

    # Backtrack to get the alignment
    alignmentS = ""
    alignmentT = ""
    i = m
    j = n

    while i > 0 and j > 0:
        if X[i][j] == X[i - 1][j - 1] + get_score(matchScore, mismatchScore, S[i - 1], T[j - 1]):
            alignmentS += S[i - 1]
            alignmentT += T[j - 1]
            i -= 1
            j -= 1
        
        elif X[i][j] == X[i - 1][j] + (cs if i > 1 and S[i - 1] == S[i - 2] else cn):
            alignmentS += S[i - 1]
            alignmentT += '-'
            i -= 1
        
        else:
            alignmentS += '-'
            alignmentT += T[j - 1]
            j -= 1

    # Add remaining characters (if any)
    while i > 0:
        alignmentS += S[i - 1]
        alignmentT += '-'
        i -= 1

    while j > 0:
        alignmentS += '-'
        alignmentT += T[j - 1]
        j -= 1

    # Reverse the aligned sequences
    alignmentS = alignmentS[::-1]
    alignmentT = alignmentT[::-1]
    
    return (X[m][n], alignmentS, alignmentT)

# Helper function to read the sequences from the fastafile using the biopython library
def read_fasta_file(filename):
    sequences = []
    for record in SeqIO.parse(filename, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def main(fasta_file_name):
    sequences = read_fasta_file(fasta_file_name)
    
    if len(sequences) < 2:
        print("Missing 1 or more sequences")
        return
    
    S = sequences[0]
    T = sequences[1]
    
    matchScore = 1
    mismatchScore = -1
    cs = -1
    cn = -2
    
    score, aligmentS, alignmentT = get_slippage_aware_alignment(S, T, matchScore, mismatchScore, cs, cn)
    
    print(f"Optimal alignment score: {score}")
    print(f"S:\n{aligmentS}")
    print(f"T:\n{alignmentT}")

if (__name__ == "__main__"):
    main("fastafile1.fa")
    main ("fastafile2.fa")
