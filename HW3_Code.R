#######################
# Homework 3 Code
#######################

#############
# Problem 3
#############

# implementing importance sampling with beta mixture model 
X = seq(0, 1, length.out = 1e2)
plot(X, dbeta(X, 5, 2))
plot(X, dbeta(X, 2, 8))

# trying to implement the mixture 
b1 <- dbeta(X, 5, 2)
b2 <- dbeta(X, 2, 8)
mix_beta = 0.3*b1 + 0.7*b2
plot(mix_beta) # mixture model, our f(x) in this case

p <- dgamma(X, 2, 5) # picking another
plot(density(mix_beta), col = "cyan") # f(x)
lines(density(p), col = "pink") # p(x)

# Now just have to multiply and find a q(x) function for importance sampling 
# taking the proudct 
pq <- p * mix_beta
plot(density(pq), col = "orange") # p(x)q(x)
# how to sample from this?

q <- dgamma(X, 3, 6)
plot(density(q))

# now we can try to get the ratio that we'll sample from 
sampling_ratio <- (p / q) * mix_beta
sampling_ratio[1] <- 0
plot(density(sampling_ratio), col = "green")

# now can do sampling...
# just to pull values from q 
# and then compute these importance values with the samplign ratio 
set.seed(100)
T = 10000
imp_vals <- rep(NA, T)

for (i in 1:T){
  y <- sample(1:100,1)
  imp_vals[i] <- sampling_ratio[y]
}

# now can get the mean 
mean_imp <- mean(imp_vals)

# probability that the random variable is included in the interval 0.45,0.55
# can do this the naive way 
target_prob <- sum(imp_vals > .45 & imp_vals < .55) / length(imp_vals)
print(paste0("Our target probability is: ", round(target_prob, 2)))
#############
# Problem 4 
#############

# Want to investigate the posterior to approximate from 
# product of a beta and log-normal random variables 
# Try simluating this in R 

mu <- dbeta(X, 2, 2)
plot(mu, type = 'l')
sigma <- dlnorm(X, 1, 10)
plot(sigma, type = 'l')

# entering the data 
Y <- c(2.3656 , 2.4952 , 1.0837 , 0.7586 , 
       0.8780 , 1.2765 , 1.4598 , 0.1801 , 
       -1.0009 , 1.4870 , -0.1193 , 0.2578)

# trying the product 
fx <- mu * sigma
plot(fx, type = 'l') # how to approximate this...
# envelope candidate function 
env <- dbeta(X, 2, 6)

# propsing this more formally 
prop_env <- function(n, p1, p2){
  return(rbeta(n, p1, p2))
}

metropolis_MCMC <- function(startval, iterations){
  chain = array(dim = c(iterations+1, 12)) # establishing chain
  chain[1,] = startvalue # length 12 along w data
  for (i in 1:iterations){
    proposal = prop_env(chain[i,], 2, 6)
    
    prob = exp(prop_env(12, 2, 6) - prop_env(chain[i,],2, 6))
    if (runif(1) < prob){
      chain[i+1, ] = proposal
    }else{
      chain[i+1, ] = chain[i,]
    }
  }
  return(chain)
}

start = Y
chain = metropolis_MCMC(startval = start , 
                        iterations = 100000)
