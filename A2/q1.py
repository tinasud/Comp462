import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# global variables
vcf_file_path = "C:\\Users\\tinas\\Downloads\\initial_population.vcf"

# Q1 a)
def read_vcf(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    data = [line.strip().split('\t') for line in lines if not line.startswith("#")]
    snp_df = pd.DataFrame(data)
    return snp_df

def initialize_population(vcf_file):
    snp_df = read_vcf(vcf_file)
    num_individuals = snp_df.shape[1] - 9 # get number of individuals

    # make a list of individuals represented by their maternal and paternal chromosomes
    population = []
    for i in range(num_individuals):
        maternal = []
        paternal = []
        for _, row in snp_df.iterrows():
            alleles = row[i + 9].split('|')  # get the alleles for individual i
            maternal.append(int(alleles[0]))
            paternal.append(int(alleles[1]))
        population.append((maternal, paternal))

    return population

def select_parent(population, fitness):
    total_fitness = sum(fitness)
    probabilities = [f / total_fitness for f in fitness]
    parent = np.random.choice(len(population), p=probabilities)
    return parent

def reproduce(parent1, parent2):
    m1, p1 = parent1
    m2, p2 = parent2
    L = len(m1)

    # pick crossover positions
    k = np.random.randint(0, L)
    k_prime = np.random.randint(0, L)

    # make child's maternal chromosome
    if np.random.rand() < 0.5:
        child_maternal = m1[:k] + p1[k:]
    else:
        child_maternal = p1[:k] + m1[k:]

    # make child's paternal chromosome
    if np.random.rand() < 0.5:
        child_paternal = m2[:k_prime] + p2[k_prime:]
    else:
        child_paternal = p2[:k_prime] + m2[k_prime:]

    return (child_maternal, child_paternal)

def simulate_generation(population, fitness):
    N = len(population)
    new_population = []
    
    for _ in range(N):
        # pick two parents
        p1_index = select_parent(population, fitness)
        p1 = population[p1_index]

        # pick second parent making sure it's different from the first
        p2_index = select_parent(population, fitness)
        while p2_index == p1_index:
            p2_index = select_parent(population, fitness)
        p2 = population[p2_index]

        # reproduce to create child
        child = reproduce(p1, p2)
        new_population.append(child)

    return new_population

def simulate_population(vcf_file, num_generations, fitness):
    population = initialize_population(vcf_file)
    
    # make sure fitness has the same length as the population
    if len(fitness) != len(population):
        fitness = [1.0] * len(population)  # reset to neutral fitness if sizes do not match

    for _ in range(num_generations):
        population = simulate_generation(population, fitness)
    return population

# to run part a uncomment the following two lines:
#initial_fitness = [1.0] * 100  # using assumption that there is a neutral fitness of 1 for each individual
#final_population = simulate_population(vcf_file_path, 10, initial_fitness)

# Q1 c)

def count_extinct_alleles(population, snp_count):
    extinct_count = 0
    for snp_index in range(snp_count):
        allele_present = False
        for individual in population:
            maternal, paternal = individual
            # check if the alternate allele is present
            if maternal[snp_index] == 1 or paternal[snp_index] == 1:
                allele_present = True
                break
        if not allele_present:
            extinct_count += 1
    return extinct_count

def estimate_extinction_probability(vcf_file, snp_count):
    population = initialize_population(vcf_file)
    population = simulate_generation(population, [1.0] * len(population))
    extinct_count = count_extinct_alleles(population, snp_count)
    extinction_prob = extinct_count / snp_count
    return extinction_prob

# to run part c uncomment the following 3 lines:
#initial_snp_count = 10000
#extinction_probability = estimate_extinction_probability(vcf_file_path, initial_snp_count)
#print(f"Estimated probability of extinction after one generation: {extinction_probability}")

# Q1 d)

def get_allele_frequencies(population, snp_count):
    frequencies = []
    for snp_index in range(snp_count):
        count = 0
        for individual in population:
            maternal, paternal = individual
            count += maternal[snp_index] + paternal[snp_index]
        frequencies.append(count / (2 * len(population)))
    return frequencies

def simulate_and_track_frequencies(vcf_file, num_generations, snp_count):
    population = initialize_population(vcf_file)
    frequency_over_time = []
    for _ in range(num_generations):
        # get frequencies before the next generation, then simulate the next generation
        frequencies = get_allele_frequencies(population, snp_count)
        frequency_over_time.append(frequencies)
        population = simulate_generation(population, [1.0] * len(population))
    return np.array(frequency_over_time)

# to run part d uncomment the following lines:
#num_generations = 20
#snp_count_to_track = 100
#allele_frequencies = simulate_and_track_frequencies(vcf_file_path, num_generations, snp_count_to_track)

#plt.figure(figsize=(12, 6))
#for i in range(snp_count_to_track):
#    plt.plot(range(num_generations), allele_frequencies[:, i], label=f'{i+1}', alpha=0.5)

#plt.title('Alternate Allele Frequencies of the First 100 SNPs Over 20 Generations')
#plt.xlabel('Generation')
#plt.ylabel('Alternate allele frequency')
#plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
#plt.show()

# Q1 e)

def get_extinction_and_fixation(population, snp_count):
    extinct_count = 0
    fixation_count = 0
    for snp_index in range(snp_count):
        allele_present = False
        homozygous = True
        for individual in population:
            maternal, paternal = individual
            if maternal[snp_index] == 1 or paternal[snp_index] == 1:
                allele_present = True
                if maternal[snp_index] == 0 and paternal[snp_index] == 0:
                    homozygous = False
            else:
                homozygous = False
        if not allele_present:
            extinct_count += 1
        elif homozygous:
            fixation_count += 1
    return extinct_count, fixation_count

def simulate_extinction_and_fixation(vcf_file, num_generations, snp_count):
    population = initialize_population(vcf_file)
    extinction_probabilities = []
    fixation_probabilities = []

    for _ in range(num_generations):
        extinct_count, fixation_count = get_extinction_and_fixation(population, snp_count)
        extinction_probabilities.append(extinct_count / snp_count)
        fixation_probabilities.append(fixation_count / snp_count)
        population = simulate_generation(population, [1.0] * len(population))

    # add 0 for both probabilities to indicate the starting generation is generation 0
    extinction_probabilities.insert(0, 0)
    fixation_probabilities.insert(0, 0)

    return extinction_probabilities, fixation_probabilities

# to run part e uncomment the following lines:
#num_generations = 1000
#snp_count = 10000
#extinction_probs, fixation_probs = simulate_extinction_and_fixation(vcf_file_path, num_generations, snp_count)

#plt.figure(figsize=(12, 6))
#plt.plot(range(num_generations + 1), extinction_probs, label='Extinction Probability', color='red')
#plt.plot(range(num_generations + 1), fixation_probs, label='Fixation Probability', color='green')
#plt.title('Probabilities of Allele Extinction and Fixation Over 1000 Generations')
#plt.xlabel('Generation')
#plt.ylabel('Probability')
#plt.legend()
#plt.show()

# Q1 f)

def calculate_fitness(individual):
    maternal, paternal = individual
    snp42_index = 42  # SNP42 index

    if maternal[snp42_index] == 1 and paternal[snp42_index] == 1:
        return 2.0  # homozygous for alternate allele
    elif maternal[snp42_index] == 1 or paternal[snp42_index] == 1:
        return 1.5  # heterozygous
    else:
        return 1.0  # all other genotypes

def select_parent_f(population):
    fitness = [calculate_fitness(ind) for ind in population]
    total_fitness = sum(fitness)
    probabilities = [f / total_fitness for f in fitness]
    parent_index = np.random.choice(len(population), p=probabilities)
    return parent_index

def simulate_generation_f(population):
    new_population = []
    N = len(population)

    for _ in range(N):
        p1_index = select_parent_f(population)
        p2_index = select_parent_f(population)
        while p2_index == p1_index:
            p2_index = select_parent_f(population)

        child = reproduce(population[p1_index], population[p2_index])
        new_population.append(child)

    return new_population

def run_simulation_f(vcf_file, num_generations, num_simulations):
    extinction_count = 0
    initial_population = initialize_population(vcf_file)

    for _ in range(num_simulations):
        population = initial_population.copy()

        for _ in range(num_generations):
            population = simulate_generation_f(population)

        # check for extinction of the alternate allele at SNP42
        snp_index = 42
        if all(individual[0][snp_index] == 0 and individual[1][snp_index] == 0 for individual in population):
            extinction_count += 1

    return extinction_count / num_simulations

# to run part f uncomment the following lines:
#num_generations = 1
#num_simulations = 1000
#extinction_probability = run_simulation_f(vcf_file_path, num_generations, num_simulations)
#print(f"Estimated probability of extinction after {num_generations} generation(s): {extinction_probability:.2f}")

# Q1 g)

def run_simulation_g(vcf_file, num_generations, num_simulations):
    extinction_count = 0
    fixation_count = 0
    initial_population = initialize_population(vcf_file)

    for _ in range(num_simulations):
        population = initial_population.copy()

        for _ in range(num_generations):
            population = [reproduce(population[select_parent_f(population)], population[select_parent_f(population)]) for _ in range(len(population))]

        # look for extinction and fixation of the alternate allele at SNP42
        snp_index = 42
        alleles = [individual[0][snp_index] for individual in population] + [individual[1][snp_index] for individual in population]
        if all(a == 0 for a in alleles):
            extinction_count += 1
        elif all(a == 1 for a in alleles):
            fixation_count += 1

    extinction_probability = extinction_count / num_simulations
    fixation_probability = fixation_count / num_simulations

    return extinction_probability, fixation_probability

# to run part g uncomment the following lines
#num_generations = 100
#num_simulations = 1000
#extinction_probability, fixation_probability = run_simulation_g(vcf_file_path, num_generations, num_simulations)
#print(f"Estimated probability of extinction after {num_generations} generations: {extinction_probability:.2f}")
#print(f"Estimated probability of fixation after {num_generations} generations: {fixation_probability:.2f}")

# Q1 h)

def calculate_fitness_h(individual):
    maternal, paternal = individual
    snp42_index = 42

    if maternal[snp42_index] == 1 and paternal[snp42_index] == 1:
        return 0.8  # homozygous for alternate deleterious allele
    elif maternal[snp42_index] == 1 or paternal[snp42_index] == 1:
        return 0.9  # heterozygous
    else:
        return 1.0  # all other genotypes

def select_parent_h(population):
    fitness = [calculate_fitness_h(ind) for ind in population]
    total_fitness = sum(fitness)
    probabilities = [f / total_fitness for f in fitness]
    parent_index = np.random.choice(len(population), p=probabilities)
    return parent_index

def run_simulation_h(vcf_file, num_generations, num_simulations):
    extinction_count = 0
    fixation_count = 0
    initial_population = initialize_population(vcf_file)

    for _ in range(num_simulations):
        population = initial_population.copy()

        for _ in range(num_generations):
            population = [reproduce(population[select_parent_h(population)], population[select_parent_h(population)]) for _ in range(len(population))]

        # look for extinction and fixation of the alternate allele at SNP42
        snp_index = 42
        alleles = [individual[0][snp_index] for individual in population] + [individual[1][snp_index] for individual in population]
        if all(a == 0 for a in alleles):
            extinction_count += 1
        elif all(a == 1 for a in alleles):
            fixation_count += 1

    extinction_probability = extinction_count / num_simulations
    fixation_probability = fixation_count / num_simulations

    return extinction_probability, fixation_probability

# to run part h uncomment the following lines:
#num_generations = 100
#num_simulations = 1000
#extinction_probability, fixation_probability = run_simulation_h(vcf_file_path, num_generations, num_simulations)
#print(f"Estimated probability of extinction after {num_generations} generations: {extinction_probability:.2f}")
#print(f"Estimated probability of fixation after {num_generations} generations: {fixation_probability:.2f}")
