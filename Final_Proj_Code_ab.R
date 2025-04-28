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

gibbs_sampler <- function(data , 
                          mean , 
                          n_iter = 5000, 
                          burnin = 1000 , 
                          global_scale = 1 , 
                          slab_scale = 1, 
                          slab_df = 1) {
  # dimensions 
  n_samples <- nrow(data)
  n_dim <- ncol(data)
  
  # storing our eventual samples 
  prec_samples <- array(0, dim = c(n_iter - burnin , n_dim, n_dim))
  
  # initializing params 
  omega <- diag(n_dim) 
  
  lambda <- matrix(1, n_dim, n_dim)
  
  # shrinkage parameter 
  tau <- global_scale
  
  # regularization component 
  c <- slab_scale^2
  
  # sum of outer products for the data 
  S <- matrix(0, n_dim, n_dim)
  for (i in 1:n_samples){
    x <- matrix(data[i, ], ncol=1)
    # fixing type 
    x <- as.numeric(x)
    S <- S + x %*% t(x)
  }
  
  # Gibbs sampling iterations 
  for (i in 1:n_iter){
    
    # regularized shrinkage coefficients 
    kappa <- c / (c + tau^2 * lambda^2)
    
    # regularized precision prior 
    prior_scale <- diag(n_dim) 
    for (j in 1:n_dim){
      for (k in 1:n_dim){
        if (j != k) # off-diagonal elements 
          prior_scale[j,k] <- kappa[j,k] * tau^2 * lambda[j,k]^2
        
      }
    }
  }
  # posterior information 
  post_df <- n_dim + n_samples 
  post_scale_inv <- solve(prior_scale) + S
  post_scale <- solve(post_scale_inv)
  
  # sampling from the wishart 
  omega <- rWishart(1, post_df, post_scale)[,,1]
  
  # updating the shrinkage parameters 
  for (j in 1:n_dim){
    for (k in 1:n_dim){
      if (j != k) {
        lambda_shape <- 1
        lambda_rate <- 1/omega[j,k]^2 + 1 / (tau^2)
        lambda[j,k] <- sqrt(1/rgamma(1, lambda_shape, lambda_rate))
        
      }
    }
  }
  
  # updating global shrinkage parameter 
  tau_shape <- (n_dim^2 - n_dim)/2 + 1
  tau_rate <- sum(1/lambda[lower.tri(lambda) | upper.tri(lambda)]^2/2 + 1/global_scale^2)
  tau <- sqrt(1/rgamma(1, tau_shape , tau_rate))
  
  # updating the slab 
  c <- 1/rgamma(1, slab_df/2, slab_df*slab_scale^2/2)
  
  if (i > burnin){
    prec_samples[i - burnin,,] <- omega 
  }
  
  # printing progress 
  if (i %% 1000 == 0){
    cat("Completed iteration", i, "of", n_iter, "\n")
  }
  
  return(list(
    precision_samples = prec_samples,
    lambda_final = lambda,
    tau_final = tau ,
    c_final = c 
  ))
}

# now testing it 
set.seed(1870)

n_dim <- 10
n_samples <- 50 

true_prec <- diag(rep(1, n_dim))

set.seed(42)
for (i in 1:5){
  j <- sample(1:n_dim, 1)
  k <- sample((1:n_dim)[-j], 1)
  val <- runif(1, 0.3, 0.8) * sample(c(-1,1),1)
  true_prec[j,k] <- true_prec[k,j] <- val
}

# ensuring positive definite 
true_prec <- true_prec + diag(n_dim) * 0.1
true_cov <- solve(true_prec)

#true_prec <- true_prec + diag(ncol(true_prec))*0.01

# now integrating our multivariate normal data 
true_mean <- rep(0, n_dim)
# we have our data 
library(MASS)
data <- mvrnorm(n_samples , true_mean,
                true_cov) # for comparison

result <- gibbs_sampler(data ,
                        true_mean,
                        n_iter = 5000,
                        burnin = 1000,
                        global_scale = 0.01,
                        slab_scale = 1,
                        slab_df = 3)

# now can try this with the proteins data...
# need the mean to be length 51
# 
n_dim <- 51
n_samples <- 256

true_mean <- rep(0, n_dim)

result_p <- gibbs_sampler(proteins.c, 
                          true_mean,
                          n_iter = 5000,
                          burnin = 1000, 
                          global_scale = 0.01, 
                          slab_scale = 1,
                          slab_df = 3)


# calculating posterior mean estimate 
precision_estimate <- apply(result_p$precision_samples, 
                            c(2,3), mean)

# thresholding for additional sparsiy as 
# mentioned in question 5 
threshold <- 0.1
precision_estimate_sparse <- precision_estimate
precision_estimate_sparse[abs(precision_estimate_sparse) < threshold] <- 0

# printing out results 
cat("True precision matrix:\n")
print(round(true_prec,2))

cat("\nEstimated precision matrix (with horseshoe prior):\n")
print(round(precision_estimate,2))

cat("\nThresholded precision matrix estimate:\n")
print(round(precision_estimate_sparse,2))

###########################
# Question 5 
# trying to make the network graph 
###########################

library(igraph)
library(ggraph)
library(tidygraph)

# function for creating a network graph 
create_dependency_graph <- function(precision_matrix, 
                                    threshold = 0, 
                                    node_names = NULL,
                                    layout = "fr") {
  # dimensions
  n_dim <- nrow(precision_matrix)
  
  # create default node names if applicable 
  if (is.null(node_names)) {
    node_names <- paste0("V", 1:n_dim)
  }
  
  # create adjacency matrix based on non-zero 
  # precision values 
  adj_mat <- matrix(0, nrow = n_dim, ncol = n_dim)
  for (i in 1:n_dim){
    for (j in 1:n_dim){
      if (i != j && abs(precision_matrix[i,j]) > threshold){
        adj_mat[i,j] <- 1
      }
    }
  }
  
  # now creating the graph 
  g <- graph_from_adjacency_matrix(adj_mat,
                                   mode = "undirected",
                                   weighted = TRUE,
                                   diag = FALSE)
  
  # setting node names 
  V(g)$name <- node_names 
  
  # edge weights on precision matrix values 
  # weights are strength of conditional dependence
  edge_list <- as_edgelist(g)
  for (i in 1:nrow(edge_list)){
    from_idx <- as.numeric(edge_list[i,1])
    to_idx <- as.numeric(edge_list[i,2])
    E(g)[from_idx %--% to_idx]$weight <- abs(precision_matrix[from_idx, to_idx])
  }
  
  return(g)
}

# function to plot the dependency graph 
plot_dependency_graph <- function(graph, 
                                  title = "Conditional Dependency Network",
                                  weight_scale = 2) {
  tg <- as_tbl_graph(graph)
  
  # creating the plot
  ggraph(tg, layout = "fr") + 
    geom_edge_link(aes(width = weight, alpha = weight),
                   color = "gray50") + 
    geom_node_point(size = 8 , color = "skyblue") + 
    geom_node_text(aes(label = name), repel = FALSE) + 
    scale_edge_width(range = c(0.5, weight_scale)) + 
    scale_edge_alpha(range = c(0.3, 0.8)) + 
    labs(title = title) + 
    theme_graph() + 
    theme(legend.position = "none")
}

# creating variable names 
var_names <- paste0("X", 1:ncol(cov_mat))

# dependency graph 
dep_graph <- create_dependency_graph(
  precision_matrix = cov_mat , 
  threshold = 0.01, 
  node_names = var_names
)

plot_dependency_graph(g)

ggsave("Conditional_density_network.png", 
       width = 10 , height = 10)

# trying to generate a labeled graph 
create_labeled_graph <- function(precision_matrix, 
                                 threshold = 0.01, 
                                 node_names = NULL,
                                 edge_labels = TRUE,
                                 title = "Conditional Dependency Network") {
  n_dim <- nrow(precision_matrix)
  
  # creating default names 
  if (is.null(node_names)){
    node_names <- paste0("V", 1:n_dim)
  }
  
  # edge list from the precision matrix 
  edges <- data.frame(from = integer(), 
                      to = integer(),
                      weight = numeric(),
                      sign = character())
  
  for (i in 1:(n_dim-1)) {
    for (j in (i+1):n_dim){
      if (abs(precision_matrix[i,j]) > threshold) {
        edges <- rbind(edges, data.frame(
          from = i, 
          to = j, 
          weight = abs(precision_matrix[i,j]),
          sign = ifelse(precision_matrix[i,j] > 0, "positive", "negative")
        ))
      }
    }
  }
  # creating nodes df 
  nodes <- data.frame(
    id = 1:n_dim , 
    label = node_names
  )
  
  # creating actual graph
  if (nrow(edges) > 0) {
    g <- tbl_graph(nodes = nodes , edges = edges)
  } else{
    g <- tbl_graph(nodes = nodes , edges = edges)
    message("Warning: No edges found above the threshold")
  }
  
  # plotting 
  p <- ggraph(g, layout = "fr") + 
    geom_edge_link(aes(width = weight , color = sign , alpha = weight),
                   show.legend = TRUE) + 
    scale_edge_color_manual(values = c("positive" = "blue", 
                                       "negative" = "red")) + 
    scale_edge_width(range = c(0.5, 3), guide = "none") + 
    scale_edge_alpha(range = c(0.5, 1), guide = "none") + 
    
    geom_node_point(size = 10 , color = "lightblue", fill = "white", shape = 21, 
                    stroke = 1.5) + 
    # node labels 
    geom_node_text(aes(label = label), repel = TRUE , size = 4,
                   fontface = "bold", bg.color = "white", bg.r = 0.15) + 
    
    # theming
    labs(title = title, 
         subtitle = paste("Edges shown: connection strength > ", threshold),
         color = "Relationship") + 
    theme_graph(base_family = "sans") + 
    theme(
      plot.title = element_text(size = 16 , face = "bold",
                                hjust = 0.5), 
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom"
    )
  
  return(p)
}

var_names <- paste0("X", 1:ncol(cov_mat))

# creating the graph
dependency_plot <- create_labeled_graph(precision_matrix =  cov_mat,
                                        threshold = 0.01,
                                        node_names = var_names , 
                                        edge_labels = TRUE)

# displaying
print(dependency_plot)


# saving 
ggsave("Conditional_density_network_labeled.png", 
       width = 10 , height = 10)
