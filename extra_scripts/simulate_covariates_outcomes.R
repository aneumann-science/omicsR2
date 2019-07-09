### Code to simulate two confounders and outcome
# Load simulated DNA methylation
data("cpg_simulated")

# Convert methylation in matrix format to data.frame
cpg_simulated.data <- as.data.frame(cpg_simulated)

# Simulate a covariate, which has a .2 standardized association with 5 CpG sites to represent a confounder
set.seed(20190704)
covariate1 <- with(cpg_simulated.data, .2*scale(cg25196575) + .2*scale(cg19821527) + .2*scale(cg01981733) + .2*scale(cg03899510) + .2*scale(cg09796913) + rnorm(500, 0, sqrt(1 - (.2^2+.2^2+.2^2+.2^2+.2^2))))

# Test, whether associations are as expected
summary(lm(covariate1 ~ scale(cg25196575) + scale(cg19821527) + scale(cg01981733) + scale(cg03899510) + scale(cg09796913), data = cpg_simulated.data))

# Simulate a second covariate, which has a .2 standardized association with 5 CpG sites to represent a confounder
set.seed(20190705)
covariate2 <- with(cpg_simulated.data, .2*scale(cg25196575) + .2*scale(cg19821527) + .2*scale(cg01981733) + .2*scale(cg03899510) + .2*scale(cg09796913) + rnorm(500, 0, sqrt(1 - (.2^2+.2^2+.2^2+.2^2+.2^2))))

# Test, whether associations are as expected
summary(lm(covariate2 ~ scale(cg25196575) + scale(cg19821527) + scale(cg01981733) + scale(cg03899510) + scale(cg09796913), data = cpg_simulated.data))

# Simulate a covariate, which has a .3 standardized association with 5 CpG sites and two covariates
set.seed(20190706)
outcome <- with(cpg_simulated.data, .3*scale(cg25196575) + .3*scale(cg19821527) + .3*scale(cg01981733) + .3*scale(cg03899510) + .3*scale(cg09796913) + .3*covariate1 + .3*covariate2 + rnorm(500, 0, sqrt(1 - (.3^2+.3^2+.3^2+.3^2+.3^2+.3^2+.3^2))))

# Test, whether associations are as expected in a model with CpG sites only
summary(lm(outcome ~ scale(cg25196575) + scale(cg19821527) + scale(cg01981733) + scale(cg03899510) + scale(cg09796913), data = cpg_simulated.data))

# Test, whether associations are as expected in a model with CpG sites and confounders
summary(lm(outcome ~ scale(cg25196575) + scale(cg19821527) + scale(cg01981733) + scale(cg03899510) + scale(cg09796913) + covariate1 + covariate2, data = cpg_simulated.data))

# Create a data.frame with the confounders and outcome
phenotype <- data.frame(ID = row.names(cpg_simulated.data), covariate1, covariate2, outcome)

# Check correlation between covariate and confounder
cor(phenotype[c("covariate1","covariate2","outcome")])

# Save data.frame with two covariates and outcome
save(phenotype, file = "data/phenotype.rda")
