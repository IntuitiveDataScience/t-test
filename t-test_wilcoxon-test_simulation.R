# Why You Should (Almost) Never Use the T-Test
# Intuitive Data Science
# MIT License

# Author(s):     Michael Crockett, Brandon Duderstadt
# Last Modified: 2019-03-09

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

set.seed(13806)

EPSILONS        <- seq(from = 0, to = 0.20, by = 0.025)
N_SAMPLES       <- 2.5*10^2
REPLICATIONS_MU <- 1.0*10^3
REPLICATIONS_SD <- 1.0*10^2
ALPHA           <- 0.05

CONTROL_SD <- 1
OUTLIER_SD <- 5
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

# Sample from distribution and get binary rejection
get_rejection <- function(sampler, n, alpha) {
  
  # Run hypothesis tests
  t_test <- t.test(sampler(n))
  w_test <- wilcox.test(sampler(n))
  
  # Store rejection outputs in vector
  rejections <- list(t = 0, w = 0)
  
  # Check if tests reject the null
  if (t_test$p.value < alpha) {
    rejections$t <- 1
  } else {
    rejections$t <- 0
  }
  
  if (w_test$p.value < alpha) {
    rejections$w <- 1
  } else {
    rejections$w <- 0
  }
  
  # Return results
  return(rejections)
  
}

###############################################################################
# PLOT DISTRIBUTION VS EPSILON
###############################################################################

# Initialize matrix to store samples
mat <- matrix(0, ncol = 2, nrow = N_SAMPLES*length(EPSILONS),
              dimnames = list(NULL, c("epsilon", "x")))

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

bp <- ggplot(df, aes(x = epsilon, y = x)) +
  geom_boxplot(colour = "black", fill = "#56B4E9") +
  scale_y_continuous(name = "x", breaks = seq(-30, 30, 10), limits = c(-30, 30)) +
  scale_x_discrete(name = "epsilon") +
  ggtitle("Samples from Control Group vs Epsilon") +
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

ggsave("synthetic-control-data.png", plot = bp, width = 6, height = 4.5, dpi = 320)

###############################################################################
# RUN HYPOTHESIS TESTS
###############################################################################

# Initialize storage vectors for mean and standard deviation
t_size  <- list(mu = rep(0, length(EPSILONS)), sd = rep(0, length(EPSILONS)))
w_size  <- list(mu = rep(0, length(EPSILONS)), sd = rep(0, length(EPSILONS)))

t_power <- list(mu = rep(0, length(EPSILONS)), sd = rep(0, length(EPSILONS)))
w_power <- list(mu = rep(0, length(EPSILONS)), sd = rep(0, length(EPSILONS)))

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(EPSILONS), style = 3)

# Loop over all epsilons
for (i in 1:length(EPSILONS)) {
  
  # Store all rejections
  t_results <- list(size_mu = rep(0, REPLICATIONS_SD), power_mu = rep(0, REPLICATIONS_SD))
  w_results <- list(size_mu = rep(0, REPLICATIONS_SD), power_mu = rep(0, REPLICATIONS_SD))
  
  # Build samplers
  center_sampler  <- sampler_factory(mu = 0, epsilon = EPSILONS[i])
  shifted_sampler <- sampler_factory(mu = SHIFTED_MU, epsilon = EPSILONS[i])
  
  # Bootstrap Monte Carlo simulation of hypothesis tests
  for (j in 1:REPLICATIONS_SD) {
    
    # Store rejections
    same_rejections <- list(t = 0, w = 0)
    diff_rejections <- list(t = 0, w = 0)
    
    # Monte Carlo simulation of hypothesis tests
    for (k in 1:REPLICATIONS_MU) {
      
      # Get rejections for same means
      rejections <- get_rejection(center_sampler, N_SAMPLES, ALPHA)
      
      same_rejections$t <- same_rejections$t + rejections$t
      same_rejections$w <- same_rejections$w + rejections$w
      
      # Get rejections for different means
      rejections <- get_rejection(shifted_sampler, N_SAMPLES, ALPHA)
      
      diff_rejections$t <- diff_rejections$t + rejections$t
      diff_rejections$w <- diff_rejections$w + rejections$w
      
    }
    
    # Store estimated size
    t_results$size_mu[j] <- same_rejections$t/REPLICATIONS_MU
    w_results$size_mu[j] <- same_rejections$w/REPLICATIONS_MU
    
    # Store estimated power
    t_results$power_mu[j] <- diff_rejections$t/REPLICATIONS_MU
    w_results$power_mu[j] <- diff_rejections$w/REPLICATIONS_MU
    
  }
  
  # Compute mean and variance for size and power
  t_size$mu[i] <- mean(t_results$size_mu)
  w_size$mu[i] <- mean(w_results$size_mu)
  t_size$sd[i] <- sd(t_results$size_mu)
  w_size$sd[i] <- sd(w_results$size_mu)
  
  t_power$mu[i] <- mean(t_results$power_mu)
  w_power$mu[i] <- mean(w_results$power_mu)
  t_power$sd[i] <- sd(t_results$power_mu)
  w_power$sd[i] <- sd(w_results$power_mu)
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
}

# Close progress bar
close(pb)

# Store results
results <- list(t_test = list(size  = list(mean = t_size$mu,
                                           sd   = t_size$sd),
                              power = list(mean = t_power$mu,
                                           sd   = t_power$sd)),
                w_test = list(size  = list(mean = w_size$mu,
                                           sd   = w_size$sd),
                              power = list(mean = w_power$mu,
                                           sd   = w_power$sd)))

file_name <- paste("monte-carlo-results_N-", N_SAMPLES,
                   "_Rmu-", REPLICATIONS_MU,
                   "_Rsd-", REPLICATIONS_SD,
                   "_Alpa-", ALPHA,
                   "_Osd-", OUTLIER_SD,
                   "_Smu-", SHIFTED_MU,
                   ".rds", sep = "")

saveRDS(results, file = file_name)

###############################################################################
# PLOT TYPE I ERROR RATE VS EPSILON
###############################################################################

df_size <- data.frame(epsilon = EPSILONS,
                      size_mu = c(results$t_test$size$mean,
                                  results$w_test$size$mean),
                      size_sd = c(results$t_test$size$sd,
                                  results$w_test$size$sd),
                      group   = c(rep("t", length(EPSILONS)),
                                  rep("w", length(EPSILONS))))

dtp <- ggplot(df_size, aes(x = as.factor(epsilon), y = size_mu, colour = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = size_mu - 2*size_sd, 
                    ymax = size_mu + 2*size_sd),
                width = 0.1) + 
  geom_line(aes(y = size_mu)) +
  scale_colour_manual(values = c("black", "red"),
                      breaks = c("w", "t"),
                      labels = c("Wilcoxon Test", "Student's T-Test")) +
  scale_y_continuous(name = "Type I Error Rate", breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_discrete(name = "epsilon") + 
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
        legend.position  = c(0.8, 0.5))

ggsave("type-i-error-vs-epsilon.png", plot = dtp, width = 6, height = 4.5, dpi = 320)

###############################################################################
# PLOT POWER VS EPSILON
###############################################################################

df_power <- data.frame(epsilon  = EPSILONS,
                       power_mu = c(results$t_test$power$mean,
                                    results$w_test$power$mean),
                       power_sd = c(results$t_test$power$sd,
                                    results$w_test$power$sd),
                       group    = c(rep("t", length(EPSILONS)),
                                    rep("w", length(EPSILONS))))

dpp <- ggplot(df_power, aes(x = as.factor(epsilon), y = power_mu, colour = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = power_mu - 2*power_sd, 
                    ymax = power_mu + 2*power_sd),
                width = 0.1) + 
  geom_line(aes(y = power_mu)) +
  scale_colour_manual(values = c("black", "red"),
                      breaks = c("w", "t"),
                      labels = c("Wilcoxon Test", "Student's T-Test")) +
  scale_y_continuous(name = "estimated power", breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_discrete(name = "epsilon") + 
  ggtitle("Estimated Power vs Epsilon") +
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
        legend.position  = c(0.8, 0.2))

ggsave("power-vs-epsilon.png", plot = dpp, width = 6, height = 4.5, dpi = 320)

###############################################################################
# CREATE ROBUSTNESS ANIMATION
###############################################################################


###############################################################################
# REAL WORLD DATASET EXAMPLE
###############################################################################

# https://www.kaggle.com/mustafaali96/weight-height

wh.data <- read.csv("/Users/mike/Documents/weight-height.csv")
wh.data$Gender <- as.factor(wh.data$Gender)

whp <- ggplot(wh.data, aes(x = Gender, y = Height)) +
  geom_boxplot(colour = "black", fill = "#56B4E9") +
  scale_y_continuous(name = "Height in Inches", breaks = seq(50, 80, 10), limits = c(50, 80)) +
  scale_x_discrete(name = "Gender") +
  ggtitle("Height by Gender") +
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

whp


