#######################
# Homework 4 Code
#######################

#######################
# Problem 1
#######################

# comparing a gibbs sampler for the centered and un-centered algos 
## loading the data 
y <- rnorm(25, 1, 1)
mean.y <- mean(y)
var.y <- var(y)
n <- length(y)

# above serve as likelihoods. Now to the gibbs sampling 
# starting values 
S <- 1000
phi <- matrix(nrow=S, ncol=2)
phi[1,] <- c(mean.y, 1/var.y)
lambda <- matrix(nrow=S, ncol=2) # for the centered case 

# Gibbs sampling algo 
set.seed(1870)
for (s in 2:S){
  # generating the new mean value 
  i <- 1
  mu_n <- rnorm(1, phi[i,1], phi[i,2] / S) # in the matrix 
  i <- i + 1
  phi[i+1,1] <- mu_n 
  
  # other parameter 
  i <- 1
  alpha_n <- rnorm(1, phi[i+1,1], phi[i,2]) # other parameter 
  phi[i+1,2] <- alpha_n
}

# we can now run this and proceed with the algorithm, 
# now we will handle the centered case 
# working with our new parameter matrix 
lambda <- matrix(nrow=S, ncol=2) # for the centered case 

lambda[1,] <- c(mean.y, mean.y + 1/var.y)

for (s in 2:S){
  # generating the new mean as before 
  i <- 1
  mu_n <- rnorm(1, lambda[i,1], lambda[i,2] / S)
  i <- i + 1
  lambda[i+1, 1] <- mu_n
  
  # other parameter 
  i <- 1
  eta_n <- rnorm(1, lambda[i+1, 1], lambda[i,2])
  lambda[i+1,2] <- eta_n 
}

#######################
# Problem 2
#######################

# reading in the glucose data 
setwd('/home/adbucks/Documents/sta_630/Homeworks/HW4/')
data <- read.table('glucose.dat')
colnames(data) <- "y"

data$y <- as.numeric(data$y)
# now can make our histogram or KDE of the data 
kde <- density(data$y)
plot(kde)

# now we can try and implement the gibbs sampler with some 
# pre-defined constraints 
a_0 <- 1
b_0 <- 1
mu_0 <- 120
tau_0 <- 200
s_0 <- 1000
rho_0 <- 10
G <- 10000 # number of iterations 

# we can start with the data 
mean_x <- mean(data$y)
var_x <- var(data$y)
n_0 <- length(data$y)

# we start with an initial beta value and then envoke our 
# mixture model 
set.seed(42)

# setting up a parameter matrix 
phi_4 <- matrix(nrow=G, ncol=3)
phi_4[1,] <- c(mean_x, 1/var_x, rho_0)


for (g in 1:G){
  # generating initial values 
  i <- 1 
  beta_s <- rbeta(1, a_0, b_0)
  i <- i + 1
  p <- beta_s
  
  # updating mu 
  mu_n <- rnorm(1, mu_0, s_0)
  mu2_n <- rnorm(1, mu_0, s_0)
  
  # not sure how to incorporate the beta 
  # and the actual mixture 
  
}

