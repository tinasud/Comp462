import pandas as pd
from scipy.stats import chi2_contingency

# global variables
genotype_file = "C:\\Users\\tinas\\Downloads\\gwas_population.vcf"
phenotype_file = "C:\\Users\\tinas\\Downloads\\gwas_phenotypes.txt"

# Q2 a)

# load genotype file
with open(genotype_file) as file:
    lines = file.readlines()
genotypes = []
for line in lines:
    if not line.startswith("#"):
        genotypes.append(line.strip().split("\t")[6:])
genotypes = pd.DataFrame(genotypes).T

# load phenotype file
phenotypes = pd.read_csv(phenotype_file, sep="\t", header=None, names=["Individual", "Disease"])

# 0 will represent homozygous reference (0|0), 1 will represent heterozygous (0|1 or 1|0), 2 will represent homozygous alternate (1|1)
genotype_df = genotypes.map(lambda x: sum(map(int, x.split('|'))))
phenotype_df = phenotypes.reset_index(drop=True)
genotype_df = genotype_df.reset_index(drop=True)

# combine along columns
combined_df = pd.concat([phenotype_df, genotype_df], axis=1).drop(['Individual'], axis=1)
contingency_df = pd.DataFrame()
p_values = []
snps = genotype_df.columns

for snp in snps:
    contingency_df = pd.crosstab(combined_df['Disease'], combined_df[snp]) # make contingency table for the current SNP
    chi2, p_value, deg_of_freedom, expected_freq = chi2_contingency(contingency_df) # perform chi-squared test of independence
    p_values.append(p_value)

results_df = pd.DataFrame({'SNP': snps, 'p_value': p_values})

# Q2 b)
print(results_df[results_df['p_value'] < 0.05])

# Q2 c)
num_of_tests = len(snps)
results_df['corrected_p_value'] = results_df['p_value'] * num_of_tests # use bonferroni correction
significant_snps_after_correction = results_df[results_df['corrected_p_value'] < 0.05]

def calculate_odds_ratio(contingency_table):
    # get counts from the contingency table
    homo_ref = contingency_table[0]
    hetero = contingency_table[1]
    homo_alt = contingency_table[2]

    # calculate probabilities for disease given each genotype
    probability_disease_homo_ref = homo_ref[1] / sum(homo_ref)
    probability_disease_hetero = hetero[1] / sum(hetero)
    probability_disease_homo_alt = homo_alt[1] / sum(homo_alt)

    odds_ratio_het = probability_disease_hetero / probability_disease_homo_ref
    odds_ratio_homo_alt = probability_disease_homo_alt / probability_disease_homo_ref
    return odds_ratio_het, odds_ratio_homo_alt

odds_ratios_of_sig_snps = []
for snp in significant_snps_after_correction['SNP']:
    contingency_table = pd.crosstab(combined_df['Disease'], combined_df[snp])
    odds_ratio_het, odds_ratio_homo_alt = calculate_odds_ratio(contingency_table)
    odds_ratios_of_sig_snps.append([snp, odds_ratio_het, odds_ratio_homo_alt])

odds_ratios_df = pd.DataFrame(odds_ratios_of_sig_snps, columns=['SNP', 'Odds_Ratio_Het', 'Odds_Ratio_Homo_Alt'])

# combine with the odd ratio df with significant snps df
final_results = pd.merge(significant_snps_after_correction[['SNP', 'corrected_p_value']], odds_ratios_df, on='SNP')
final_results = pd.merge(results_df[['SNP', 'p_value']], final_results, on='SNP')
print(final_results)
