#################################
# Final Project Code 
#################################

# importing data 
setwd('/home/adbucks/Documents/sta_630/Final_Project/')

# reading from excel 
library(readxl)
df <- read_xls('aml-rppa.xls')

#################################
# Question 1 
# Normalize the data 
# to ensure the assumption mu = 0 holds
#################################

# want to start by isolating the 51 proteins 
library(dplyr)

# now grabbing just the proteins, everything after col 47
proteins <- df[ ,48:98] # checked the index 

head(proteins) # we know these are all normally distributed 
# trying this another way 
proteins.c <- apply(proteins, MARGIN = 2 , function(x) scale(x, scale = FALSE))

# converting 
proteins.c <- as.data.frame(proteins.c)
# can now center to ensure that mu = 0 
# ensuring they're all basically 0 
colMeans(proteins.c) < .05 # checks out 

#################################
# Question 2
# Based on the joint regression model , 
# derive the full conditional distributions 
# for each off-diagonal element ij and each diagonal
# element of the precision matrix 
# specify appropriate prior distributions as 
# part of your derivation 
#################################

# might want to start by getting the covariance matrix and 
# by extension the precision matrix , can interpret these in the following way 
cov_mat <- cov(proteins.c)
prec_mat <- 1 / cov(proteins.c)

# seems as if this model framework wants to combine 
# two pseudo-likelihoods and picking a prior 
# for a given precision matrix entry , done to preserve symmetry in the 
# precision matrix 
# given that we are actually estimating the precision matrix 
# can derive the full conditional from the assignment, along with picking an 
# appropriate prior for the precision matrix 

# given Gelman and other relevant literature, we will pick the Wishart as a 
# conjugate prior for the precision matrix when deriving the full conditional 

#################################
# Question 3

# reflect the sparsity assumption of the 
# underlying graph by incorporating
# sparsity-inducing priors, like those given in the 
# hw 
# pick a prior and justify your choice 
#################################

# deciding between the spike and slab prior or the adaptive shrinkage 
# prior 

# after reviewing the literature (can add this to references) will be using 
# and implementing a regularized horseshoe prior, given thaqt the spike and 
# slab is sensitive to slot width and this provides a greater ease of 
# computation 

# now just implementing this into the full conditional 
# we will stick with the same model framework, and will introduce the prior
# using the half cauchy distribution 

# we can say the beta values in the linear model have the normal prior with 
# mean 0 and the variance signified by tau squared and lambda squared, with a 
# regularization term defined for the lambda parameter 
# the prior for the lambda value is the half-cauchy

# $\beta_j | \lambda_j , \tau , c \sim N(0, \tau^2, \lambda^2)$
# \lambda_j \sim C^+(0,1)
# \lambda_j^2 = \frac{c^2 \lambda_j^2}{c^2 + \tau^2 \lambda_j^2}

#################################
# Question 4

# Implement a Gibbs sampler 
# based on the earlier full conditional
# distributions 
#################################

# want to think about how the full conditionals might work with each other
# before implementing the gibbs sampler 

  

  
