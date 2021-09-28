
This code was used in [this preprint](https://doi.org/10.1101/856898).

# Dependencies
- numpy
- pandas
- matplotlib
- seaborn
- scipy

summary=../../human_gwas_data/summary_statistics/${f%_summary_statistics.txt}_summary_statistics.txt
out=../../human_gwas_data/results_2021/${f%_summary_statistics.txt}_model
python wcrep.py $summary $out


# Step 1: Training models

Both the winner's curse model and the winner's curse and confounding models are trained using the wcrep.py python script

## Usage
python wcrep.py [data table] [prefix]

## Input data
The input data is as follows:
1. A table with the z-scores (tab deliminated with a header; the columns can be anything)
	- The first column should have z-scores of the significant variants in the discovery study
	- The second column should have the corresponding z-scores in the replication study
	- The third column should have the sample size of the first study
	- The fourth column should have the size of the second study
	- The threshold used for the discovery study
2. Output prefix

## Output data
The output is a single text file with the following values in each column.
1. Name of the input data
2. Sample size of the discovery study
3. Sample size of the replication study
4. Discovery threshold
1. Observed replication rate
2. Expected replication rate under the winner's curse only model
3. Expected replication rate under the winner's curse and confounding model
4. Standard deviation parameter for the true genetic effect size
5. Standard deviation parameter for the study-specific confounding in the first study
6. Standard deviation parameter for the study-specific confounding in the second study
7. Number of variants significant in the discovery study



