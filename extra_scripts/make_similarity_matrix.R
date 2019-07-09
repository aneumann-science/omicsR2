# Load BGData for similarity matrix computation
library(BGData)

# Load methylation data
data("cpg_simulated")

# Compute similarity matrix for DNA methyaltion
# getG does the following:
# 1. z-score standardizes the methylation values
# 2. Tatakes the product of the transpose
# 3. Scales the similarity matrix so diagonal is 1
Gmt <- getG(cpg_simulated)

# Save methylation similarity matrix
save(Gmt, file = "data/Gmt.rda")

