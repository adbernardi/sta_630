#################################
# Final Project Code v2
#################################

# importing data 
setwd('/home/adbucks/Documents/sta_630/Final_Project')

library(readxl)
df <- read_xls('aml-rppa.xls') # reads as a tibble 

# trimming the data to get only the proteins
library(dplyr)

proteins <- df[, 48:98] # just getting the MVN proteins 
# checking distribution 
plot(density(proteins$ACTB)) # basically right 

# now centering so that the mean can be 0 
proteins.c <- apply(proteins, MARGIN = 2 , function(x) scale(x, scale = FALSE))

proteins.c <- as.data.frame(proteins.c) # type conversion 

# checking 
sum(apply(proteins.c, MARGIN = 2 , mean) < 0.05) == ncol(proteins.c) # true 

#################################
# full conditionals and the gibbs sampler 
#################################

# setting conditions for the algorithm
n <- nrow(proteins.c)
p <- ncol(proteins.c)

# establishing results from the data/evidence
Sigma_obs <- cov(proteins.c)

# more initializing 
Omega <- diag(p)
Sigma <- solve(Omega)

# trying a storage matrix 
Omega_samples <- array(0, dim = c(p, p, 5000)) # iterations minus burn in 

# running the gibbs sampler with the conditionals 
# given for the high-dimension MVN 

# loading relevant libs 
library(mvtnorm)
library(Matrix)

# setting hyperparams 
lambda <- 1
# trying a function over a loop 
gibbs <- function(data,
                  n_iter = 10000, 
                  burn_in = 5000,
                  p, 
                  n,
                  Sigma
){
  for (iter in 1:n_iter){
    for (i in 1:p){
      not_i <- (1:p)[-1]
      
      S_11 <- Sigma[not_i, not_i]
      S_12 <- Sigma[not_i, i]
      
      Omega_11 <- Omega[not_i , not_i]
      beta_mean <- -solve(Omega_11) %*% S_12
      beta_cov <- solve(Omega_11) / n
      
      beta <- t(rmvnorm(1, mean = as.numeric(beta_mean),
                        sigma = beta_cov))
      omega_ii <- rgamma(1, shape = (n / 2) + 1, 
                         rate = (Sigma_obs[i,i] + t(beta) %*% Omega_11 %*% beta) / 2)
      Omega[i,i] <- omega_ii # filling in 
      Omega[not_i, i] <- Omega[i,not_i] <- as.numeric(beta) # flattening
    }
    # lasso style shrinkage 
    Omega[-diag(p)] <- Omega[-diag(p)] / (1 + lambda)
    
    # storing 
    if (iter > burn_in){
      Omega_samples[,,iter - burn_in] <- Omega
    }
    
    # progress bar
    if (iter %% 1000 == 0) cat("Iteration", iter, "\n")
  }
  # actually returning stuff 
  #return(Omega_post_mean)
  return(Omega_samples)
}

# calling 
set.seed(1870)

Omega_est <- gibbs(proteins.c, p = p , n = n , Sigma = Sigma_obs)

# can now do inference and move on to graphing from here
Omega_post_mean <- apply(Omega_est , c(1,2), mean) # not actually sure why
# the mean might not make sense here 

# getting partial correlations as per the assignment 
partial_corr <- -Omega_post_mean / sqrt(pmax(0, diag(Omega_post_mean) %o% diag(Omega_post_mean)))

diag(partial_corr) <- 1

# now can get into thresholding and graphing 
# can change the threshold for differing effects 
threshold <- 0.1
adjacency_mat <- (abs(partial_corr) > threshold) * 1
diag(adjacency_mat) <- 0

# maybe a quick fix 
adjacency_mat <- ifelse(adjacency_mat != 0 && adjacency_mat != 1, 0, adjacency_mat)

library(igraph)
# creating the actual graph object 
g <- graph_from_adjacency_matrix(adjacency_mat , 
                                 mode = "Undirected",
                                 weighted = TRUE,
                                 diag = FALSE)

# labels created, now just ned to maybe make label font size
# smaller and spread the nodes out somewhat 
plot(g, vertex.label = colnames(proteins.c)) #seems to work!!

library(ggplot2)
# saving plot 


#################################
# MCMC diagnostics and further analysis 
#################################

library(coda)
library(bayesplot)
# can use gelman methods here?
# a type of trace plot 
# type conversion for diagnositcs 
# have to drop the NAs before proceeding 
Omega_est[is.na(Omega_est)] <- 0

mcmc_trace(Omega_est)
# need to change parameter names, so probably need to do
# something with changing the column names to be in the 
# same order as the actual proteins, like for the 
# undirected graph



