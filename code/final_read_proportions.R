


################################################################################
# Evaluating the proportion of endophytes the raw reads got. To determine 
# the abundance. Compare the estimated vs measured
################################################################################


props <- read.table("raw_reads_prop.tsv", sep = "\t",
                      header = TRUE)
samples <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
             "X08", "X090", "X095")
props$samples <- samples


props[,-6] * 100

hc <- c(0, 0.25, 0.43, 0.57, 0.66, 0.76, 0.81, 0.88, 0.92, 0.96, 0.98)

ratio <- c(0, 2.539, 2.146, 1.888, 1.660, 1.519, 1.363, 1.255, 1.154, 1.071, 1.034)

OG_hc <- c(seq(0,0.9, by = 0.1), 0.95)



cor.test(OG_hc, ratio, method = "pearson")


# Perform linear regression
model <- lm(ratio ~ hc)

# Extract coefficients
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Print the equation
cat("y =", slope, "* x +", intercept)
# y = -0.06264557 * x + 1.461936




# Perform linear regression
model <- lm(OG_hc ~ ratio)

# Extract coefficients
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Print the equation
cat("y =", slope, "* x +", intercept)
# y = -0.1133251 * x + 0.6564689
