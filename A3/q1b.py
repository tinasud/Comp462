import pandas as pd
import numpy as np
from Bio import SeqIO
from math import log

def viterbi(sequence, config_file): 
    neg_inf = float("-inf")
    states = ["inter", "start", "middle", "stop"]
    sequence_length = len(sequence)
    
    # initialize probability and pointer matrices
    probabilities = np.full((4, sequence_length), float("-inf"))
    pointer = np.full((4, sequence_length), -9)
    
    emissions_intergenic = []
    emissions_start = []
    emissions_middle = []
    emissions_stop = []
    
    emissions_intergenic_df = pd.DataFrame()
    emissions_start_df = pd.DataFrame()
    emissions_middle_df = pd.DataFrame()
    emissions_stop_df = pd.DataFrame()
    
    # read the config file to parse the emission and transition probabilities
    with open(config_file, "r") as c:
        lines = c.readlines()
        
        # identify the indices of lines that correspond to each state (start, stop, middle)
        start_index, stop_index, middle_index = 0, 0, 0
        for i, line in enumerate(lines):
            if "start" in line: 
                start_index = i
            elif "stop" in line: 
                stop_index = i
            elif "middle" in line: 
                middle_index = i
        
        # read emissions for the intergenic state
        for i in range(11, 15):
            temp = lines[i].split(' ')
            # avoid zero emissions by replacing them with a very small value (0.000001)
            emissions_intergenic.append(float(temp[1].strip('\n')) if float(temp[1].strip('\n')) != 0 else 0.000001)
        # store intergenic emissions in a dataFrame, indexed by nucleotide (A, T, G, C)
        emissions_intergenic_df = pd.DataFrame(emissions_intergenic, index=["A", "T", "G", "C"], columns=["intergenic emission"])
        
        emissions_start = [(line.split(' ')[0].strip(':'), float(line.split(' ')[1].strip('\n'))) 
                           for line in lines[start_index + 1: stop_index]]
        emissions_start_df = pd.DataFrame(emissions_start, columns=["start codon", "start emission"]).set_index("start codon")
        
        emissions_stop = [(line.split(' ')[0].strip(':'), float(line.split(' ')[1].strip('\n'))) 
                          for line in lines[stop_index + 1: middle_index]]
        emissions_stop_df = pd.DataFrame(emissions_stop, columns=["stop codon", "stop emission"]).set_index("stop codon")
        
        emissions_middle = [(line.split(' ')[0].strip(':'), float(line.split(' ')[1].strip('\n'))) 
                            for line in lines[middle_index + 1:]]
        emissions_middle_df = pd.DataFrame(emissions_middle, columns=["middle codon", "middle emission"]).set_index("middle codon")
        
        avg_intergenic_length = float(lines[3].split(' ')[5].strip('\n'))
        avg_cds_length = float(lines[8].split(' ')[5].strip('\n'))
    
    # define transition probabilities between states
    transmissions = pd.DataFrame(
        {"inter": [((avg_intergenic_length - 1) / avg_intergenic_length), 
                   (1 / avg_intergenic_length), 0, 0],  # inter -> inter, inter -> start, inter -> middle, inter -> stop
        "start": [0.0, 0.0, 1.0, 0.0],  # start -> start, start -> middle
        "middle": [0.0, 0.0, ((avg_cds_length - 1) / avg_cds_length), (1 / avg_cds_length)],  # middle -> middle, middle -> stop
        "stop": [1.0, 0.0, 0.0, 0.0]},  # stop -> intergenic
        index=states)

    # set the initial probabilities for the first nucleotide
    probabilities[0, 0] = float(emissions_intergenic_df.loc[sequence[0]].iloc[0])  # set intergenic state probability at position 0
    pointer[0, 0] = 0  # point to intergenic state
    
    # fill viterbi table
    for i in range(1, sequence_length):
        # compute probabilities for the intergenic state at position i
        temp_list = [
            probabilities[0, i - 1] + log(transmissions.loc["inter", "inter"]),  # intergenic -> intergenic transition
            probabilities[3, i - 3] + log(transmissions.loc["inter", "stop"])  # stop -> intergenic transition
        ]
        probabilities[0, i] = log(emissions_intergenic_df.loc[sequence[i]].values[0]) + max(temp_list)  # add emission probability
        
        # store pointer to the state with the highest probability
        if temp_list[0] > temp_list[1] and pointer[0, i - 1] == 0:
            pointer[0, i] = 0  # point to intergenic state
        if temp_list[0] < temp_list[1] and i >= 3:
            if pointer[3, i - 3] == 2 and max(probabilities[:, i - 3]) == probabilities[3, i - 3]:
                pointer[0, i] = 3  # point to stop state if conditions are met
        
        # handle transitions for the start, stop, and middle states
        if i <= sequence_length - 3:
            triplets = sequence[i] + sequence[i + 1] + sequence[i + 2] # get the current triplet (3 nucleotides)
            
            # stop codon emissions and transition
            if triplets in emissions_stop_df.index and probabilities[2, i - 3] != neg_inf:
                probabilities[3, i] = log(emissions_stop_df.loc[triplets].values.item()) + probabilities[2, i - 3] + log(transmissions.loc["stop", "middle"])
                if pointer[2, i - 2] != 1 and probabilities[2, i - 1] != 1:
                    pointer[3, i] = 2  # point to middle state if conditions are met
            
            # start codon emissions and transition
            elif triplets in emissions_start_df.index and probabilities[0, i - 1] != neg_inf:
                probabilities[1, i] = log(emissions_start_df.loc[triplets].values.item()) + probabilities[0, i - 1] + log(transmissions.loc["start", "inter"])
                if max(probabilities[:, i - 1]) == probabilities[0, i - 1]:
                    pointer[1, i] = 0  # point to intergenic state if conditions are met

            # middle codon emissions and transition
            if triplets in emissions_middle_df.index and (probabilities[1, i - 3] != neg_inf or probabilities[2, i - 3] != neg_inf):
                temp_list = [
                    probabilities[1, i - 3] + log(transmissions.loc["middle", "start"]),  # start -> middle transition
                    probabilities[2, i - 3] + log(transmissions.loc["middle", "middle"])  # middle -> middle transition
                ]
                probabilities[2, i] = log(emissions_middle_df.loc[triplets].iloc[0].item()) + max(temp_list)  # add middle emission probability
                if temp_list[0] > temp_list[1] and pointer[1, i - 3] == 0:
                    pointer[2, i] = 1  # point to start state
                elif temp_list[0] < temp_list[1] and pointer[2, i - 3] in [1, 2]:
                    pointer[2, i] = 2  # point to middle state

    return probabilities, pointer

def backtrack(sequence_list, probabilities, pointer):
    max_prob = float('-inf')
    curr_index = 0
    start_indices = []
    stop_indices = []

    # find the state of the last nucleotide with the highest probability
    for i in range(4):  # loop through the possible 4 states (0, 1, 2, 3)
        if probabilities[i, len(sequence_list) - 1] > max_prob:
            # update the max probability and set the current state to the one with the highest probability
            max_prob = pointer[i, (len(sequence_list) - 1)]
            curr_index = i

    # start backtracking from the last nucleotide until we reach the begininng of the sequence
    x = len(sequence_list) - 1
    while x >= 0:
        if curr_index == 0:  # intergenic : state 0
            # check if we can transition from state 3 at position (x-3) to state 0 at position x
            if (pointer[3, x - 3] == 2) and (pointer[0, x] == 3):  # transition conditions for state 0
                stop_indices.append(x)  # mark this position as a stop
                curr_index = 3  # update state to 3
                x -= 3  # skip the next 2 nucleotides
            else:
                x -= 1  # otherwise, move back 1 nucleotide
                
        elif curr_index == 2:  # middle : state 2
            # check for a specific transition condition from state 1 to state 2
            if pointer[1, x - 3] == 0 and (pointer[2, x - 3] == -9):  # condition to move from start to middle
                start_indices.append(x - 2)  # mark the start position of this segment
                curr_index = 1  # transition to state 1
                x = x - 3  # skip back 3 nucleotides
            else:
                x -= 3  # otherwise, move back 3 nucleotides
                
        elif curr_index == 1:  # state 1: another state (likely corresponding to a region transition)
            # check if we need to transition to intergenic (state 0)
            if pointer[2, x - 3] == -9:  # condition to switch to state 0 (intergenic)
                curr_index = 0  # update to intergenic state
                x -= 1  # move back by 1 nucleotide
            else:
                curr_index = 2  # otherwise, remain in state 2
                x -= 3  # skip back 3 nucleotides
        
        else:  # default case: if current state is neither 0, 1, nor 2, set state to 2 and skip back 3
            curr_index = 2  # set to state 2 (middle state)
            x = x - 3  # skip back 3 nucleotides
    
    start_indices.reverse()
    stop_indices.reverse()
    return start_indices, stop_indices

def write_gff3_file(contig, start, stop):
    for i in range (0, len(start), 1):
        with open("results.gff3","a+") as f:
            newline = contig + "\tena\tCDS\t" + str(start[i]) + "\t" + str(stop[i]) + "\t.\t+\t0\t.\t\n"
            f.writelines(newline)

def main(config_file, fasta_file):
    fasta_list = list(SeqIO.parse(fasta_file, "fasta"))
    
    for i in range (0, len(fasta_list), 1):
        probabilities, pointer = viterbi(fasta_list[i].seq, config_file)
        start, stop = backtrack(list(fasta_list[i].seq), probabilities, pointer)
        write_gff3_file(fasta_list[i].id, start, stop)

main("config.txt", "C:\\Users\\tinas\\Downloads\\Vibrio_vulnificus.ASM74310v1.dna.nonchromosomal.fa")
