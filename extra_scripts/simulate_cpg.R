# Load winsorized DNA methylation at birth
load("GENR_cord_meth_3IQR_Winsorized.RData")

# Select 1000 CpG sites randomly
birth_win_5000 <- as.data.frame(birth_win[,sample(1:ncol(birth_win),1000)])

# Extract mean and SD for the 1000 CpG sites
# Simulate a random distribution with the mean and SD of the actual measurement
# n= 500 participants
# Add some random noise (max/min +- 10% of SD)
cpg_simulated <- sapply(birth_win_5000, function(probe) {
	probe_mean <- mean(probe)
	probe_sd <- sd(probe)
	rnorm(500, probe_mean, probe_sd)+runif(500,-probe_sd/10,probe_sd/10)
})

# Simulate participant ID
row.names(cpg_simulated) <- paste0(sample(1:10000, 500), "_R01C01")

# Save the simulated DNA methylation dataset
save(cpg_simulated, file = "cpg_simulated.rda")
