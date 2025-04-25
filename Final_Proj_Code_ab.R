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

# can now center to ensure that mu = 0 
proteins.c <- sapply(proteins, function(x) scale(x, scale = FALSE))

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

