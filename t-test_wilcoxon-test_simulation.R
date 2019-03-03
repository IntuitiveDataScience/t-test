# Why You Should (Almost) Never Use the T-Test
# Intuitive Data Science

# Author(s):     Michael Crockett, Brandon Duderstadt
# Last Modified: 2019-03-03

###############################################################################
# LOAD LIBRARIES
###############################################################################

library(sysfonts)
library(xkcd)
library(ggplot2)
library(animation)

###############################################################################
# DEFINE SIMULATION PARAMETERS
###############################################################################

set.seed(20190303)

EPSILONS     <- seq(from = 0, to = 0.20, by = 0.025)
N_SAMPLES    <- 2.5*10^2
REPLICATIONS <- 1.0*10^3
ALPHA        <- 0.05

OUTLIER_SD <- 10
SHIFTED_MU <- 0.25

###############################################################################
# DEFINE FUNCTIONS
###############################################################################

# Creates sampler for mixed normal distribution
sampler_factory <- function(mu, epsilon) {
  
  # Draws n samples from mixed normal distribution
  sampler <- function(n) {
    
    # Use uniform random variables on [0, 1] and epsilon threshold
    sapply(runif(n), FUN = function(u) {
      
      # Increasing epsilon increases samples from outlier distribution
      if (u < epsilon) {
        rnorm(1, mean = mu, sd = OUTLIER_SD)
      } else {
        rnorm(1, mean = mu, sd = CONTROL_SD)
      }
      
    })
  }
  
  return(sampler)
  
}

###############################################################################
# PLOT DISTRIBUTION VS EPSILON
###############################################################################

# Initialize matrix to store samples
mat <- matrix(0, ncol = 2, nrow = N_SAMPLES*length(EPSILONS),
              dimnames = list(NULL, c("epsilon", "y")))

# Store samples drawn with specified epsilon
for (i in 1:length(EPSILONS)) {
  
  sampler <- sampler_factory(mu = 0, epsilon = EPSILONS[i])
  
  mat[i*(1:N_SAMPLES), 1] <- EPSILONS[i]
  mat[i*(1:N_SAMPLES), 2] <- sampler(N_SAMPLES)
  
}

# Convert matrix to dataframe
df         <- as.data.frame(mat)
df$epsilon <- as.factor(df$epsilon)

# Build boxplot
fill <- "#4271AE"
line <- "#1F3552"

bp <- ggplot(df, aes(x = epsilon, y = y)) +
  geom_boxplot(colour = "black", fill = "#56B4E9") +
  scale_y_continuous(name = "y", breaks = seq(-40, 40, 10), limits = c(-40, 40)) +
  scale_x_discrete(name = "epsilon") +
  ggtitle("Boxplot of Synthetic Outlier Data") +
  theme(axis.line.x      = element_line(size = 0.5, colour = "black"),
        axis.line.y      = element_line(size = 0.5, colour = "black"),
        axis.line        = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_blank(),
        panel.background = element_blank(),
        plot.title       = element_text(size = 20, family = "xkcd", hjust = 0.5),
        text             = element_text(size = 16, family = "xkcd"),
        axis.text.x      = element_text(colour = "black", size = 12),
        axis.text.y      = element_text(colour = "black", size = 12))

bp

###############################################################################
# RUN HYPOTHESIS TESTS
###############################################################################

# Store test results from Monte Carlo simulation
same_test_size  <- list(t = rep(0, length(EPSILONS)), w = rep(0, length(EPSILONS)))
same_test_power <- list(t = rep(0, length(EPSILONS)), w = rep(0, length(EPSILONS)))

diff_test_size  <- list(t = rep(0, length(EPSILONS)), w = rep(0, length(EPSILONS)))
diff_test_power <- list(t = rep(0, length(EPSILONS)), w = rep(0, length(EPSILONS)))

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(EPSILONS), style = 3)

# Loop over all epsilons
for (i in 1:length(EPSILONS)) {
  
  # Build samplers
  center_sampler  <- sampler_factory(mu = 0, epsilon = EPSILONS[i])
  shifted_sampler <- sampler_factory(mu = SHIFTED_MU, epsilon = EPSILONS[i])
  
  # Log correct and incorrect rejections
  same_incorrect_rejections    <- list(t = 0, w = 0)
  different_correct_rejections <- list(t = 0, w = 0)
  
  # Monte Carlo simulation of hypothesis tests
  for (j in 1:REPLICATIONS) {
    
    # Generate samples
    centered_w_outliers  <- center_sampler(N_SAMPLES)
    shifted_w_outliers   <- shifted_sampler(N_SAMPLES)
    
    # Run t-test
    p_t_same      <- t.test(centered_w_outliers)$p.value
    p_t_different <- t.test(shifted_w_outliers)$p.value
    
    # Run Wilcoxon test
    p_w_same      <- wilcox.test(centered_w_outliers)$p.value
    p_w_different <- wilcox.test(shifted_w_outliers)$p.value
    
    # Store t-test result for same means
    if (p_t_same < ALPHA) {
      same_incorrect_rejections$t <- same_incorrect_rejections$t + 1
    }
    
    # Store t-test result for different means
    if (p_t_different < ALPHA) {
      different_correct_rejections$t <- different_correct_rejections$t + 1
    }
    
    # Store Wilcoxon result for same means
    if (p_w_same < ALPHA) {
      same_incorrect_rejections$w <- same_incorrect_rejections$w + 1
    }
    
    # Store Wilcoxon result for different means
    if (p_w_different < ALPHA) {
      different_correct_rejections$w <- different_correct_rejections$w + 1
    }
    
  }
  
  # Store test size for particular epsilon
  same_test_size$t[i] <- same_incorrect_rejections$t/REPLICATIONS
  same_test_size$w[i] <- same_incorrect_rejections$w/REPLICATIONS
  
  # Store test power for particular epsilon
  diff_test_power$t[i] <- different_correct_rejections$t/REPLICATIONS
  diff_test_power$w[i] <- different_correct_rejections$w/REPLICATIONS
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
}

# Close progress bar
close(pb)

###############################################################################
# PLOT TYPE I ERROR RATE VS EPSILON
###############################################################################

df_size <- data.frame(epsilon     = EPSILONS,
                      same_t_size = same_test_size$t,
                      same_w_size = same_test_size$w)

dtp <- ggplot(df_size, aes(x = epsilon)) +
  geom_line(aes(y = same_t_size, colour = "Student's t-test")) +
  geom_point(aes(y = same_t_size), colour = "black") +
  geom_line(aes(y = same_w_size, colour = "Wilcoxon Test")) +
  geom_point(aes(y = same_w_size), colour = "red") +
  scale_colour_manual(values = c("black", "red"), 
                      breaks = c("Wilcoxon Test", "Student's t-test")) +
  scale_y_continuous(name = "type i error rate", breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_continuous(name = "epsilon", breaks = EPSILONS, limits = range(EPSILONS)) +
  ggtitle("Type I Error Rate vs Epsilon") +
  theme(axis.line.x      = element_line(size = 0.5, colour = "black"),
        axis.line.y      = element_line(size = 0.5, colour = "black"),
        axis.line        = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_blank(),
        panel.background = element_blank(),
        plot.title       = element_text(size = 20, family = "xkcd", hjust = 0.5),
        text             = element_text(size = 16, family = "xkcd"),
        axis.text.x      = element_text(colour = "black", size = 12),
        axis.text.y      = element_text(colour = "black", size = 12),
        legend.title     = element_blank(),
        legend.text      = element_text(size = 16, family = "xkcd"),
        legend.position  = c(0.8, 0.55))

dtp

###############################################################################
# PLOT POWER VS EPSILON
###############################################################################

df_power <- data.frame(epsilon      = EPSILONS,
                       diff_t_power = diff_test_power$t,
                       diff_w_power = diff_test_power$w)

dpp <- ggplot(df_power, aes(x = epsilon)) +
  geom_line(aes(y = diff_t_power, colour = "Student's t-test")) +
  geom_point(aes(y = diff_t_power), colour = "black") +
  geom_line(aes(y = diff_w_power, colour = "Wilcoxon Test")) +
  geom_point(aes(y = diff_w_power), colour = "red") +
  scale_colour_manual(values = c("black", "red"), 
                      breaks = c("Wilcoxon Test", "Student's t-test")) +
  scale_y_continuous(name = "power", breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_continuous(name = "epsilon", breaks = EPSILONS, limits = range(EPSILONS)) +
  ggtitle("Power vs Epsilon") +
  theme(axis.line.x      = element_line(size = 0.5, colour = "black"),
        axis.line.y      = element_line(size = 0.5, colour = "black"),
        axis.line        = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_blank(),
        panel.background = element_blank(),
        plot.title       = element_text(size = 20, family = "xkcd", hjust = 0.5),
        text             = element_text(size = 16, family = "xkcd"),
        axis.text.x      = element_text(colour = "black", size = 12),
        axis.text.y      = element_text(colour = "black", size = 12),
        legend.title     = element_blank(),
        legend.text      = element_text(size = 16, family = "xkcd"),
        legend.position  = c(0.8, 0.5))

dpp
