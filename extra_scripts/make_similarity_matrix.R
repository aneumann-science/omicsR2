# Load BGData for similarity matrix computation
library(BGData)

# Load methylation data
data("cpg_simulated")

# Compute similarity matrix for DNA methyaltion
# getG does the following:
# 1. z-score standardizes the methylation values
# 2. Takes the product of the transpose
# 3. Scales the similarity matrix so diagonal is 1
Gmt <- getG(cpg_simulated)

# Save methylation similarity matrix
save(Gmt, file = "data/Gmt.rda")

# Create matrix representing batch effects
# Every 500 participants are assigned to one batch
# This represents proper randomization to batches
Batch <- diag(500)
Batch[1:100,1:100] <- 1
Batch[100:200,100:200] <- 1
Batch[200:300,200:300] <- 1
Batch[300:400,300:400] <- 1
Batch[400:500,400:500] <- 1
save(Batch, file = "data/Batch.rda")

