
This code was used in [this preprint](https://doi.org/10.1101/856898).

# Dependencies
Training (R):

-mvtnorm

Visualization (python):

-numpy
-pandas
-matplotlib
-scipy

# Step 1: Training models

Both the winner's curse model and the winner's curse and confounding models are trained using the model_replication.R script.  

## Input data
The input data is as follows:
1. A table with the z-scores
	- The first column should have z-scores of the significant variants in the discovery study
	- The second column should have the corresponding z-scores in the replication study
	- The table should be tab deliminated with a header (The columns can be named anything)
2. The sample size of the first study
3. The sample size of the first study
4. A path for an output text file containing model parameters

## Output data
The output is a single text file with the following values in each line.
1. The observed replication rate
2. The expected replication rate under the winner's curse only model
3. The expected replication rate under the winner's curse and confounding model
4. The variance parameter for the true genetic effect size
5. The variance parameter for the study-specific confounding in the first study
6. The variance parameter for the study-specific confounding in the second study

## Sample usage

Rscript src/model_replication.R sample/input.txt 252972 80022 sample/output.txt 

# Step 2: Visualizing models

An example of how to visualize the distributions learned is shown in the plot_models.py script.  The input to this script is as follows.
1. The model output file learned in Step 1
2. The input data used to train the models
3. The sample size of the first study
4. The sample size of the second study
5. The path of the output plot
6. A title for the output plot

Example usage:
python src/plot_models.py sample/output.txt sample/input.txt 252972 80022 sample/output.png test


